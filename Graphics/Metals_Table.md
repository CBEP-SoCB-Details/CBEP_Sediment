Developing a Table for Metals
================
Curtis C. Bohlen, Casco Bay Estuary Partnership
October 22, 2020

  - [Introduction](#introduction)
  - [Load Libraries](#load-libraries)
  - [Load Data](#load-data)
      - [Folder References](#folder-references)
      - [Metals Data](#metals-data)
          - [Change Factor Levels](#change-factor-levels)
      - [Select Metals to Be Analyzed](#select-metals-to-be-analyzed)
      - [Estimating Non-detects](#estimating-non-detects)
  - [Robust Regression models](#robust-regression-models)
      - [T-S Slopes](#t-s-slopes)
      - [Check Significance with Kendall’s
        Tau](#check-significance-with-kendalls-tau)
  - [Generating a Summary Table of
    Metals](#generating-a-summary-table-of-metals)
      - [Trend Words](#trend-words)
      - [Region Model](#region-model)
          - [Examine Pairwise
            Comparisons](#examine-pairwise-comparisons)
          - [Hand-generated comparison
            text](#hand-generated-comparison-text)

<img
    src="https://www.cascobayestuary.org/wp-content/uploads/2014/04/logo_sm.jpg"
    style="position:absolute;top:10px;right:50px;" />

# Introduction

Casco Bay Estuary Partnership collected sediment samples in 1991, 1994,
2000, 2001, 2002, 2010, and 2011 to look at concentrations of toxic
contaminants in Casco Bay surface Sediments. These studies were
complemented by data collected by under the auspices of EPA’s the
National Coastal Assessment (NCA) and National Coastal Condition
Assessment (NCCA).

Chemicals studied included a number of different metals. Here we develop
a narrative table for the State of the Bay Report.

# Load Libraries

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.3     v dplyr   1.0.2
    ## v tidyr   1.1.2     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(mblm)
library(emmeans)

library(CBEPgraphics)
load_cbep_fonts()
theme_set(theme_cbep())

library(LCensMeans)
```

# Load Data

## Folder References

``` r
sibfldnm <- 'Derived_Data'
parent   <- dirname(getwd())
sibling  <- file.path(parent,sibfldnm)
niecefldnm <- 'Data_Subsets'
niece <- file.path(sibling,niecefldnm)
fn <- "metals.csv"
```

## Metals Data

``` r
metals_data <- read_csv(file.path(niece,fn),
                      col_types = cols(.default = col_character(),
                                        Sample_Year = col_double(),
                                        Replicate = col_integer(),
                                        CASRN = col_skip(),
                                        Result = col_double(),
                                        MDL = col_double(),
                                        RL = col_double(),
                                        Det_Flag = col_integer(),
                                        Qualifier = col_skip(),
                                        `QA Qualifier` = col_skip(),
                                        Reportable_Result = col_skip(),
                                        ERL = col_double(),
                                        ERM = col_double() )
                      ) %>%
  mutate(Replicate = Replicate == -1) %>%
  mutate(Det_Flag = Det_Flag == 1) %>%
  mutate(nd_flag = ! Det_Flag) %>%
  mutate(Parameter = if_else(Parameter == "Chromium (total)",
                             "Chromium",
                             Parameter))
```

### Change Factor Levels

``` r
metals_data <- metals_data %>%
  mutate(LVL = factor(LVL, levels = c('Below ERL','Between ERL and ERM',
                                     'Above ERM'))) %>%
  mutate(Region = factor(Region, levels = c("Inner Bay",
                                            "West Bay",
                                            "East Bay",
                                            "Outer Bay",
                                            "Cape Small"))) %>%
  mutate(Era = ordered(Era, levels = c( "1990s", "2010s", "2000s")))
```

## Select Metals to Be Analyzed

We restrict analysis to those metals sampled in all three Eras.

``` r
xt <- xtabs(~ Parameter + Era, data = metals_data)
rowcount <- apply(xt,1,function(x) sum(x>0,na.rm = TRUE))
selected <- names(rowcount[rowcount>2])

metals_data <- metals_data %>%
  filter(Parameter %in% selected)
rm(selected, xt)
```

## Estimating Non-detects

Using our lognormal maximum likelihood procedure. We apply this
procedure to the entire data set for each metal at one time, which may
bias results if data should not be considered as being drawn from a
single distribution (e.g., if concentrations are dramatically different
in different Regions or Eras). Because detection limits are uniformly
low, often well below most observed values, this is likely to have
relatively small effect on slope estimates, which are our focus here.

Note that conditional means associated with rare non-detects often can
not be readily estimated using the existing algorithm, because it
requires drawing a huge sample to get enough samples below the detection
limits.

We up the parameter `max_samp` to 5 million in hopes of reducing the
number of times that happens. That slows this computation so that it
takes several minutes, but even so, we get a few observations where we
can not estimate a conditional mean. (Note that each notification,
below, is for a single observation). This behavior may be addressed in a
future version of the LCMeans package.

This is a time-consuming step, so we set cache = true. That reduces time
for reknitting the notebook, but does not help when running the notebook
interactively. An alternative would be to save a version of the data
(with estimated conditional means) after this step, and not recalculate,
but that runs the risk of leaving in outdated data if upstream data is
modified.

``` r
metals_data <- metals_data %>%
  mutate(val = if_else(nd_flag, MDL, Result)) %>%
  group_by(Parameter) %>%
  mutate(Res_ML = sub_cmeans(val, nd_flag, max_samp = 5*10^6)) %>%
  ungroup(Parameter)
```

    ## Estimated sample size required >5e+06. Returning NA.
    ## 
    ## Estimated sample size required >5e+06. Returning NA.
    ## 
    ## Estimated sample size required >5e+06. Returning NA.
    ## 
    ## Estimated sample size required >5e+06. Returning NA.
    ## 
    ## Estimated sample size required >5e+06. Returning NA.

# Robust Regression models

We can conduct a robust regression, based on Theil-Sen slope estimators
(actually a slightly more robust estimator by Siegel).

We could not get the following to work inside a pipe or `lapply()` call,
so we fell back on using a loop. Also, mblm does not like having a log
transform in the model formula, so we had to move the log transform
outside the call.

``` r
mods <-metals_data %>%
  group_by(Parameter) %>%
  nest()

lmods <- list()
for (n in seq_along(mods$Parameter)) {
  metal <- mods$Parameter[n]
  tmp <- metals_data %>%
    filter(Parameter == metal) %>%
    mutate(logval = log(Res_ML)) %>%
    filter(! is.na(Res_ML))
  rlm <- mblm(logval ~ Sample_Year, data = tmp)
  lmods[[metal]] <- rlm
}

mods$rlm <- lmods
rm(lmods, rlm, tmp, metal)
```

## T-S Slopes

``` r
rlm_slopes <- unlist(lapply(mods$rlm ,function(m) coefficients(m)[[2]]))
names(rlm_slopes) <- mods$Parameter
```

## Check Significance with Kendall’s Tau

There is some controversy regarding the appropriateness of different
approaches to estimating statistical significance of these slopes. The
default approach provided in the `mblm` package is not well supported.
Instead, we use the related Kendall’s Tau, which is generally more
conservative.

``` r
mods <- mods %>%
  mutate(tau_test = 
           lapply(data,
                  function(df) cor.test(~ log(Res_ML) + Sample_Year,
                                        data = df,
                                        method = 'kendall')))
```

``` r
tau_p_vals = lapply(mods$tau_test, function(t) t$p.value)
names(tau_p_vals) <- mods$Parameter

sig_tau <- names(tau_p_vals[tau_p_vals<0.05])
rlm_slopes[sig_tau]
```

    ##      Cadmium     Chromium         Lead       Nickel         Zinc 
    ##  0.023733256 -0.008904760 -0.007323551 -0.008605040 -0.006956498

So, we now end up with cadmium showing an increase, Chromium, Lead, and
Nickle declining, and other metals showing no trend.

# Generating a Summary Table of Metals

A graphic for metals will take up too much room in the report, so we
want to generate a table summarizing results for metals. We want to
indicate “Increaseing”, “Decreasing”, “No Trend” for each metal, and add
something about the fraction of observations exceeding thresholds in
2010.

## Trend Words

``` r
tau =   unlist(lapply(mods$tau_test, function(t) t$estimate))
tau_words <- with(mods, if_else(! (tau_p_vals < 0.05),
                             "No Trend",
                             if_else(rlm_slopes < 0, 'Decreasing',
                                     'Increasing')))
names(tau_words) <- mods$Parameter
tau_words
```

    ##      Arsenic      Cadmium     Chromium       Copper         Iron         Lead 
    ##   "No Trend" "Increasing" "Decreasing"   "No Trend"   "No Trend" "Decreasing" 
    ##       Nickel     Selenium         Zinc      Mercury       Silver 
    ## "Decreasing"   "No Trend" "Decreasing"   "No Trend"   "No Trend"

Remember, we are ONLY working with the metals for which we have data
from all three eras here. Also, we look at the ERL because no metal
samples in 2010 exceeded the ERM values.

The categories are defined as follows: \* Never: Never observed
(frequency = 0)  
\* Rare: \< 10%  
\* Uncommon: \< 25%  
\* Common: \<50%  
\* Frequent \> 50%

``` r
tab1 <- metals_data %>%
  filter(Era == '2010s') %>%
  group_by(Parameter) %>%
  summarize(has_ERL = any(! is.na(ERL)),
            num  = sum(! is.na(Result)),
            num_exceeds = sum(Result > ERL, na.rm = TRUE),
            pct_exceeds = round(num_exceeds/num,3)*100,
            labs = cut(pct_exceeds,
                       c(-0.1, 0.1, 10,25, 50, 110),
                       labels = c('Never',
                                      'Rare',
                                      'Uncommon',
                                      'Common',
                                      'Frequent')),
            `Exceeds ERL†` = if_else(has_ERL, as.character(labs), 'No ERL')) %>%
  select(-has_ERL, -num, -num_exceeds, -pct_exceeds, -labs)
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
tab1
```

    ## # A tibble: 11 x 2
    ##    Parameter `Exceeds ERL†`
    ##    <chr>     <chr>         
    ##  1 Arsenic   Frequent      
    ##  2 Cadmium   Rare          
    ##  3 Chromium  Common        
    ##  4 Copper    Rare          
    ##  5 Iron      No ERL        
    ##  6 Lead      Rare          
    ##  7 Mercury   Common        
    ##  8 Nickel    Frequent      
    ##  9 Selenium  No ERL        
    ## 10 Silver    Never         
    ## 11 Zinc      Never

Now, we want to do something similar, adding a columns indicating the
direction of any SIGNIFICANT slope estimates to the table. We already
worked these out based on the linear regression models, above.

Here’s the code, slightly modified, to avoid a long-distance dependency
in the R Notebook that could prove confusing.

``` r
is_sign   <- unlist(lapply(mods$rlm,function(m) summary(m)$coefficients[2,4] < 0.025))
decreasing <- unlist(lapply(mods$rlm,function(m) coefficients(m)[2] < 0))

tab2 <- tibble(Parameter = mods$Parameter,
               is_sign   = is_sign,
               is_decreasing = decreasing,
               'Trend' = if_else(! is_sign, 'No Trend',
                              if_else(is_decreasing,
                                      'Decreasing','Increasing'))) %>%
  select(-is_sign, -is_decreasing)
```

And last, lets identify whether there are significant differences among
Regions

## Region Model

``` r
mods <- mods %>%
  mutate(region_anova = lapply(data,
                      function(df) lm(log(Res_ML) ~ Region, data = df)))
```

Here we pull the p value from the ANOVA table, rather than the summary,
because we want to an overall test for inequality across all regions.
’summary.lm()`provides an array, but`anova.lm()\` produces a data
frame, with non-syntactic variable names.

See `Metals_Analysis_3.Rmd` for more in-depth look at methods and
assumptions.

``` r
region_p_vals = lapply(mods$region_anova, function(t) anova(t)$`Pr(>F)`[1])
names(region_p_vals) <- mods$Parameter

(sig_region <- names(region_p_vals[region_p_vals<0.025]))
```

    ##  [1] "Arsenic"  "Cadmium"  "Chromium" "Copper"   "Iron"     "Lead"    
    ##  [7] "Nickel"   "Zinc"     "Mercury"  "Silver"

So almost all metals show differences by region.

### Examine Pairwise Comparisons

In `Metals_Analysis~3.Rmd`, we used the related `pwpp()` function for a
graphical depiction of pairwise comparisons. We skip that here to save
space.

``` r
compare <- list()
for (n in seq_along(mods$Parameter)) {
  metal <- mods$Parameter[n]
  compare[[metal]] <- emmeans(mods$region_anova[[n]],
                              ~ Region, type = 'Response')
}
mods$region_emm <- compare
rm(metal, compare)

for (n in seq_along(mods$Parameter)) {
  cat('\n---------------------', mods$Parameter[n], '---------------------\n')
  print(mods$region_emm[[n]])
}
```

    ## 
    ## --------------------- Arsenic ---------------------
    ##  Region     emmean     SE  df lower.CL upper.CL
    ##  Inner Bay    2.25 0.0626 225     2.13     2.38
    ##  West Bay     2.29 0.0647 225     2.16     2.41
    ##  East Bay     2.32 0.0752 225     2.18     2.47
    ##  Outer Bay    2.21 0.0727 225     2.07     2.35
    ##  Cape Small   1.82 0.1075 225     1.61     2.03
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## --------------------- Cadmium ---------------------
    ##  Region     emmean     SE  df lower.CL upper.CL
    ##  Inner Bay  -0.861 0.0885 224    -1.04   -0.687
    ##  West Bay   -1.112 0.0907 224    -1.29   -0.934
    ##  East Bay   -0.986 0.1054 224    -1.19   -0.778
    ##  Outer Bay  -1.691 0.1019 224    -1.89   -1.491
    ##  Cape Small -2.903 0.1508 224    -3.20   -2.606
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## --------------------- Chromium ---------------------
    ##  Region     emmean     SE  df lower.CL upper.CL
    ##  Inner Bay    4.16 0.0449 225     4.07     4.24
    ##  West Bay     4.22 0.0464 225     4.13     4.31
    ##  East Bay     4.30 0.0539 225     4.19     4.40
    ##  Outer Bay    4.17 0.0521 225     4.06     4.27
    ##  Cape Small   3.92 0.0771 225     3.77     4.07
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## --------------------- Copper ---------------------
    ##  Region     emmean     SE  df lower.CL upper.CL
    ##  Inner Bay    3.00 0.0611 225     2.88     3.12
    ##  West Bay     2.76 0.0632 225     2.64     2.89
    ##  East Bay     2.79 0.0734 225     2.64     2.93
    ##  Outer Bay    2.49 0.0709 225     2.35     2.63
    ##  Cape Small   1.73 0.1050 225     1.52     1.93
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## --------------------- Iron ---------------------
    ##  Region     emmean     SE  df lower.CL upper.CL
    ##  Inner Bay    10.2 0.0415 225     10.2     10.3
    ##  West Bay     10.3 0.0429 225     10.2     10.4
    ##  East Bay     10.3 0.0498 225     10.2     10.4
    ##  Outer Bay    10.2 0.0482 225     10.1     10.3
    ##  Cape Small   10.0 0.0713 225      9.9     10.2
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## --------------------- Lead ---------------------
    ##  Region     emmean     SE  df lower.CL upper.CL
    ##  Inner Bay    3.71 0.0338 225     3.65     3.78
    ##  West Bay     3.25 0.0350 225     3.18     3.32
    ##  East Bay     3.27 0.0406 225     3.19     3.35
    ##  Outer Bay    3.27 0.0393 225     3.20     3.35
    ##  Cape Small   2.93 0.0581 225     2.82     3.05
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## --------------------- Nickel ---------------------
    ##  Region     emmean     SE  df lower.CL upper.CL
    ##  Inner Bay    3.19 0.0480 225     3.10     3.29
    ##  West Bay     3.20 0.0497 225     3.10     3.30
    ##  East Bay     3.23 0.0577 225     3.12     3.34
    ##  Outer Bay    3.12 0.0558 225     3.01     3.23
    ##  Cape Small   2.77 0.0825 225     2.61     2.94
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## --------------------- Selenium ---------------------
    ##  Region     emmean    SE  df lower.CL upper.CL
    ##  Inner Bay  -0.837 0.155 216   -1.142   -0.532
    ##  West Bay   -1.091 0.159 216   -1.404   -0.779
    ##  East Bay   -0.581 0.184 216   -0.944   -0.217
    ##  Outer Bay  -0.971 0.180 216   -1.326   -0.616
    ##  Cape Small -1.458 0.323 216   -2.095   -0.822
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## --------------------- Zinc ---------------------
    ##  Region     emmean     SE  df lower.CL upper.CL
    ##  Inner Bay    4.45 0.0477 225     4.35     4.54
    ##  West Bay     4.29 0.0494 225     4.19     4.38
    ##  East Bay     4.34 0.0573 225     4.23     4.45
    ##  Outer Bay    4.18 0.0554 225     4.07     4.29
    ##  Cape Small   3.75 0.0820 225     3.59     3.91
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## --------------------- Mercury ---------------------
    ##  Region     emmean     SE  df lower.CL upper.CL
    ##  Inner Bay   -1.72 0.0833 217    -1.88    -1.56
    ##  West Bay    -2.88 0.0855 217    -3.05    -2.71
    ##  East Bay    -2.29 0.0993 217    -2.48    -2.09
    ##  Outer Bay   -2.78 0.0970 217    -2.97    -2.59
    ##  Cape Small  -3.74 0.1681 217    -4.07    -3.40
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## --------------------- Silver ---------------------
    ##  Region     emmean     SE  df lower.CL upper.CL
    ##  Inner Bay   -1.13 0.0827 225    -1.29   -0.966
    ##  West Bay    -1.57 0.0855 225    -1.74   -1.405
    ##  East Bay    -1.56 0.0993 225    -1.75   -1.360
    ##  Outer Bay   -1.87 0.0960 225    -2.06   -1.677
    ##  Cape Small  -2.54 0.1421 225    -2.82   -2.258
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95

### Hand-generated comparison text

We could not figure out an easy way to generate a summary description of
the pairwise comparisons, so we generated text by hand, and read them
into a data frame here so we can access them later to construct a table

``` r
tab3 <- tribble(
~Parameter     , ~`Comparison of Regions`,
'Arsenic'    , 'Cape Small is low.',
'Cadmium'    , 'Outer Bay is low;  Cape Small is even lower.',
'Chromium'   , 'Cape Small is lower than East Bay and West Bay.',
'Copper'     , 'Outer Bay is low; Cape Small even lower.',
'Iron'       , 'Cape Small is lower than East Bay and West Bay.',
'Lead'       , 'Cape Small is low;  Inner Bay is high.',
'Mercury'    , 'Cape Small < Outer Bay and West Bay < East Bay < Inner Bay',
'Nickel'     , 'Cape Small is low.',
'Selenium'   , 'No differences.',
'Silver'     , 'Cape Small is low; Inner Bay is high.',
'Zinc'       , 'Cape Small is low. Outer Bay lower than inner Bay.')
```

``` r
tab <- tab1 %>% left_join(tab2)  %>% left_join(tab3) %>%
  rename(Metal = Parameter)
```

    ## Joining, by = "Parameter"
    ## Joining, by = "Parameter"

``` r
cap = '† Rare: < 10%; Uncommon: < 25%; Common: <50%, Frequent > 50%'
knitr::kable(tab, caption = cap)
```

| Metal    | Exceeds ERL† | Trend      | Comparison of Regions                                         |
| :------- | :----------- | :--------- | :------------------------------------------------------------ |
| Arsenic  | Frequent     | Decreasing | Cape Small is low.                                            |
| Cadmium  | Rare         | Increasing | Outer Bay is low; Cape Small is even lower.                   |
| Chromium | Common       | Decreasing | Cape Small is lower than East Bay and West Bay.               |
| Copper   | Rare         | Decreasing | Outer Bay is low; Cape Small even lower.                      |
| Iron     | No ERL       | Decreasing | Cape Small is lower than East Bay and West Bay.               |
| Lead     | Rare         | Decreasing | Cape Small is low; Inner Bay is high.                         |
| Mercury  | Common       | No Trend   | Cape Small \< Outer Bay and West Bay \< East Bay \< Inner Bay |
| Nickel   | Frequent     | Decreasing | Cape Small is low.                                            |
| Selenium | No ERL       | Decreasing | No differences.                                               |
| Silver   | Never        | Increasing | Cape Small is low; Inner Bay is high.                         |
| Zinc     | Never        | Decreasing | Cape Small is low. Outer Bay lower than inner Bay.            |

† Rare: \< 10%; Uncommon: \< 25%; Common: \<50%, Frequent \> 50%
