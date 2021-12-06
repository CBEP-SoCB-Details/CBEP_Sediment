Review CBEP Historical Sediment Data
================
Curtis C. Bohlen, Casco Bay Estuary Partnership
October 13, 2020

  - [Introduction](#introduction)
  - [Load Libraries](#load-libraries)
  - [Load Data](#load-data)
      - [Folder References](#folder-references)
  - [Sums Data](#sums-data)
  - [Examine Structure of Sums Data](#examine-structure-of-sums-data)
  - [Units](#units)

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

Chemicals studied included metals, polycyclic aromatic hydrocarbons
(PAHs), polychlorinated biphenyls (PCBs), organochlorine pesticides,
dioxins and furans, and organotins. These contaminants are all
persistent in the marine environment.

# Load Libraries

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------------------------------------------------------- tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.3     v dplyr   1.0.2
    ## v tidyr   1.1.2     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts ------------------------------------------------------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(readxl)
library(knitr)

# library(GGally)

library(CBEPgraphics)
load_cbep_fonts()

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
fn <- "sums_totals.csv"
```

# Sums Data

``` r
sums_data <- read_csv(file.path(niece,fn),
                      col_types = cols(.default = col_character(),
                                        Sample_Year = col_double(),
                                        Replicate = col_integer(),
                                        CASRN = col_skip(),
                                        Result = col_double(),
                                        MDL = col_skip(),
                                        RL = col_skip(),
                                        Det_Flag = col_skip(),
                                        Qualifier = col_skip(),
                                        `QA Qualifier` = col_skip(),
                                        Reportable_Result = col_skip(),
                                        ERL = col_double(),
                                        ERM = col_double() )
                      ) %>%
  mutate(Replicate = Replicate == -1)
```

# Examine Structure of Sums Data

``` r
xtabs(~Parameter + Era, data = sums_data)
```

    ##                         Era
    ## Parameter                1990s 2000s 2010s
    ##   Butyltin (mono+di+tri)    32    32    59
    ##   CDD/CDF (total)           30    32    17
    ##   Chlordane (total)         65    76    82
    ##   DDT+DDE+DDD (sum)         65    76    82
    ##   PAHs (High MW)            65    78    82
    ##   PAHs (Low MW)             65    78    82
    ##   PAHs (total)              65    78    82
    ##   PCBs (total)              65    76    82
    ##   Pesticides (total)        65    76    82

But MANY of those observations are zeros, because of high detection
limits in the last couple of years. Those appear in these data as NAs:

``` r
xtabs(~Parameter + Era, data = sums_data, subset = ! is.na(Result))
```

    ##                         Era
    ## Parameter                1990s 2000s 2010s
    ##   Butyltin (mono+di+tri)    31    19     6
    ##   CDD/CDF (total)           30    32    17
    ##   Chlordane (total)         64    12     0
    ##   DDT+DDE+DDD (sum)         65    75     1
    ##   PAHs (High MW)            65    78    82
    ##   PAHs (Low MW)             65    78    80
    ##   PAHs (total)              65    78    82
    ##   PCBs (total)              65    76     1
    ##   Pesticides (total)        65    76     1

So non-detects are a real problem for Pesticides, PCBs, and Butyltins in
2000 and 2001.

**The metadata is not entirely clear, but it appears (from examining the
Access database) that these sums omit non-detects, effectively equating
non-detects to zero. That is inconsistent with how we handled
non-detects in several other toxics data sets, where we have been using
maximum likelihood estimators of expected values.**

# Units

``` r
xtabs(~Parameter  + Era + Units , data = sums_data)
```

    ## , , Units = ng/g dry
    ## 
    ##                         Era
    ## Parameter                1990s 2000s 2010s
    ##   Butyltin (mono+di+tri)    32    32    59
    ##   CDD/CDF (total)            0     0     0
    ##   Chlordane (total)         65    76    82
    ##   DDT+DDE+DDD (sum)         65    76    82
    ##   PAHs (High MW)            65    78    82
    ##   PAHs (Low MW)             65    78    82
    ##   PAHs (total)              65    78    82
    ##   PCBs (total)              65    76    82
    ##   Pesticides (total)        65    76    82
    ## 
    ## , , Units = ng/kg dry
    ## 
    ##                         Era
    ## Parameter                1990s 2000s 2010s
    ##   Butyltin (mono+di+tri)     0     0     0
    ##   CDD/CDF (total)           30    32    17
    ##   Chlordane (total)          0     0     0
    ##   DDT+DDE+DDD (sum)          0     0     0
    ##   PAHs (High MW)             0     0     0
    ##   PAHs (Low MW)              0     0     0
    ##   PAHs (total)               0     0     0
    ##   PCBs (total)               0     0     0
    ##   Pesticides (total)         0     0     0

Ramboll Standardized units in the Access database, so, MOST sums are
expressed in ng/g dry weight (\~ ppb).

The Dioxins and Furans are expressed in ng/kg, or pg/g or approximately
parts per trillion.

There are no Squirts for Dioxins and Furans. Instead, Ramboll ALSO
expressed them in TEQ – Tox equivalents. Toxic equivalents provide a way
to estimate the cumulative toxic effect of a mixture of related
chemicals by weighting each compound by its relative toxic effect,
compared to some reference compound (conventionally TCDD).
