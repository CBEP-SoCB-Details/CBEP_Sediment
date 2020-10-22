Review CBEP Historical Sediment Data
================
Curtis C. Bohlen, Casco Bay Estuary Partnership
October 13, 2020

  - [Introduction](#introduction)
      - [Sample Locations](#sample-locations)
  - [Load Libraries](#load-libraries)
  - [Load Data](#load-data)
      - [Folder References](#folder-references)
      - [Load Data](#load-data-1)
          - [Add the “Era” variable](#add-the-era-variable)
  - [Sums Data](#sums-data)
  - [Examine Data Structure](#examine-data-structure)
      - [Categories Sampled Each Year](#categories-sampled-each-year)
      - [Parameters and Parameter
        Groups](#parameters-and-parameter-groups)
  - [Units](#units)
  - [Subdivide Data By Parameter
    Group](#subdivide-data-by-parameter-group)
  - [Check for Other Sums and Totals in the
    Data](#check-for-other-sums-and-totals-in-the-data)
      - [The “Official” Totals](#the-official-totals)
      - [Checking All Data Subsets](#checking-all-data-subsets)
          - [Butyl Tin Totals](#butyl-tin-totals)
          - [Dioxin Totals](#dioxin-totals)
          - [Conclusion](#conclusion)
      - [Metals Data Totals](#metals-data-totals)
      - [PAH Data Totals](#pah-data-totals)
      - [PCB Data Totals](#pcb-data-totals)
      - [Pesticides Data Totals](#pesticides-data-totals)
      - [Physical Data Totals](#physical-data-totals)
          - [Checks on sums and totals](#checks-on-sums-and-totals)
      - [Total Grain Size](#total-grain-size)
  - [Examine Structure of Sums Data](#examine-structure-of-sums-data)
  - [Parameters and Parameter
    Groups](#parameters-and-parameter-groups-1)
  - [Units](#units-1)
  - [Load SQuiRTs Data](#load-squirts-data)
      - [Units](#units-2)

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

## Sample Locations

The original (1991) sampling locations were selected by an expert group
by placing dots on a map of Casco Bay to get sample locations
distributed evenly in areas known to not have rocky substrate (and
where, therefore, soft sediment samples could be collected). Originally,
each location was given a location code, of the form “XX\#\#”, where
“XX” is a two letter region code, and “\#\#” is a two digit sequence
number.

An effort was made in subsequent years to sample the same locations, but
for a variety of reasons, samples were not always collected in quite the
same locations. In the Early 2010s, CBEP staff realized sample locations
were recorded inconsistently, and made a concerted effort to identify
where each sample had actually been collected, often by going back to
paper field data sheets associated with sampling events.

A large fraction, but not all, observations were collected within a
hundred meters or so of nominal original (1991) sample locations. We
grouped sample locations into Locations (clusters of samples within
about 100m of the original sampling locations), assigning new Location
Codes as needed. Within Locations we assigned Substations (within the
resolution of the GPS equipment used at the time, typically \~ 5 to 10
meters, which is smaller than expected drift during sample collection.).

To simplify data analysis and visualization, samples are “binned” into
five major subregions of the Bay, as follows

(descriptions from 2017 Ramboll Environ report to CBEP):

> Inner Bay: The Inner Bay includes the western most part of Casco Bay,
> encompassing the Fore River and Presumpscot River estuaries. The
> cities of Portland and South Portland and the towns of Falmouth and
> Yarmouth are located within its drainage.

> Outer Bay: The Outer Bay includes the large open water area to the
> south of the other regions. It represents the area that connects Casco
> Bay to the rest of the Gulf of Maine.

> West Bay: West Bay extends from Yarmouth on the west to Orrs and
> Bailey Islands on the east. It includes Maquoit and Middle Bays,
> Harpswell Sound, and the communities of Freeport, Harpswell, and parts
> of southern Brunswick. The Royal and Harraseeket Rivers discharge into
> West Bay.

> East Bay: East Bay includes the inland portions of Casco Bay bordered
> by Orrs and Bailey Islands on the west and Phippsburg on the east. It
> includes Quahog Bay and the New Meadows River in southern Brunswick
> and Phippsburg.

> Cape Small: Cape Small is the eastern most region of the bay. It
> includes the southern end of the Phippsburg peninsula to Small Point.
> The mouth of the Lower Kennebec River flows into the Gulf of Maine to
> the east of Phippsburg. While it is not part of Casco Bay, coastal
> circulation patterns indicate that the discharge from the Lower
> Kennebec (which includes flows from the Kennebec and Androscoggin
> Rivers) is entrained into the bay in the Cape Small area (Janzen et
> al. 2005).

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
sibfldnm <- 'Original_Data'
parent   <- dirname(getwd())
sibling  <- file.path(parent,sibfldnm)
fn <- "Casco_Bay_Sediment_Data_Query.xlsx"
```

## Load Data

``` r
sedtox_data <- read_excel(file.path(sibling,fn),
    col_types = c("text", "text", "text", 
        "text", "numeric", "numeric", "text", 
        "numeric", "date", "numeric", "text", 
        "text", "text", "text", "text", "numeric", 
        "numeric", "numeric", "text", "text", 
        "numeric", "text", "text", "text", 
        "text"))
```

### Add the “Era” variable

``` r
sedtox_data <- sedtox_data %>%
  mutate(Era = (Sample_Year-1990) %/% 10)%>%
  mutate(Era = factor(Era, labels = c('1990s', '2000s', '2010s')))

xtabs(~Sample_Year + Era, data = sedtox_data)
```

    ##            Era
    ## Sample_Year 1990s 2000s 2010s
    ##        1991 10203     0     0
    ##        1994   818     0     0
    ##        2000     0  2402     0
    ##        2001     0  3816     0
    ##        2002     0  2357     0
    ##        2010     0     0  6182
    ##        2011     0     0  6594

# Sums Data

``` r
fn <- "Sums_Query.xlsx"

sums_data <- read_excel(file.path(sibling,fn), 
    col_types = c("text", "text", "text", 
        "text", "numeric", "numeric", "text", 
        "numeric", "date", "numeric", "text", 
        "text", "text", "text", "text", "numeric", 
        "numeric", "numeric", "text", "text", 
        "numeric", "numeric", "numeric", 
        "text", "text")) %>%
  mutate(Era = (Sample_Year-1990) %/% 10)%>%
  mutate(Era = factor(Era, labels = c('1990s', '2000s', '2010s')))
```

# Examine Data Structure

## Categories Sampled Each Year

``` r
xtabs(~Sample_Year + Parameter_Group, data = sedtox_data)
```

    ##            Parameter_Group
    ## Sample_Year CDDF Inorganic Organotin  PAH  PCB Pesticide Physical Radio
    ##        1991    0       715         0 2015 5200      1950      323     0
    ##        1994  690         0       128    0    0         0        0     0
    ##        2000    0       405         0  756  525       575      141     0
    ##        2001  437       825        76  924  693       759       83    19
    ##        2002  299       585        52  504  378       414      106    19
    ##        2010  459       864        85 1120 1936      1210      508     0
    ##        2011    0       966       210 1176 2268      1344      630     0

The category “Radio” includes only measurement of Uranium-238, not
otherwise studied, so of little value. We don’t track it further.

## Parameters and Parameter Groups

``` r
the_order <- as.numeric(factor(sedtox_data$Parameter_Group))

tmp <- sedtox_data %>%
  mutate(Parameter = factor(Parameter,
                            levels = levels(fct_reorder(sedtox_data$Parameter, the_order))))
                                                 
knitr::kable(xtabs(~Parameter + Parameter_Group, data = tmp))
```

|                                             | CDDF | Inorganic | Organotin | PAH | PCB | Pesticide | Physical | Radio |
| :------------------------------------------ | ---: | --------: | --------: | --: | --: | --------: | -------: | ----: |
| 1,2,3,4,6,7,8-Heptachlorodibenzo-p-dioxin   |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 1,2,3,4,6,7,8-Heptachlorodibenzofuran       |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 1,2,3,4,7,8-Hexachlorodibenzo-p-dioxin      |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 1,2,3,4,7,8-Hexachlorodibenzofuran          |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 1,2,3,4,7,8,9-Heptachlorodibenzofuran       |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 1,2,3,6,7,8-Hexachlorodibenzo-p-dioxin      |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 1,2,3,6,7,8-Hexachlorodibenzofuran          |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 1,2,3,7,8-Pentachlorodibenzo-p-dioxin       |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 1,2,3,7,8-Pentachlorodibenzofuran           |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 1,2,3,7,8,9-Hexachlorodibenzo-p-dioxin      |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 1,2,3,7,8,9-Hexachlorodibenzofuran          |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 2,3,4,6,7,8-Hexachlorodibenzofuran          |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 2,3,4,7,8-Pentachlorodibenzofuran           |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 2,3,7,8-Tetrachlorodibenzo-p-dioxin         |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| 2,3,7,8-Tetrachlorodibenzofuran             |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| CDD/CDF (total)                             |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| HpCDD (total)                               |   17 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| HpCDF (total)                               |   17 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| HxCDD (total)                               |   17 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| HxCDF (total)                               |   17 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| Octachlorodibenzo-p-dioxin                  |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| Octachlorodibenzofuran                      |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| PCB-077 (3,3’,4,4’-TeCB) (dioxin-like)      |   62 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| PCB-126 (3,3’,4,4’,5-PeCB) (dioxin-like)    |   62 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| PCB-169 (3,3’,4,4’,5,5’-HxCB) (dioxin-like) |   62 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| PeCDD (total)                               |   17 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| PeCDF (total)                               |   17 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| TCDD (total)                                |   17 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| TCDF (total)                                |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| TEQ CDD/CDF                                 |   79 |         0 |         0 |   0 |   0 |         0 |        0 |     0 |
| Aluminum                                    |    0 |       165 |         0 |   0 |   0 |         0 |        0 |     0 |
| Antimony                                    |    0 |       165 |         0 |   0 |   0 |         0 |        0 |     0 |
| Arsenic                                     |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Barium                                      |    0 |       113 |         0 |   0 |   0 |         0 |        0 |     0 |
| Beryllium                                   |    0 |       113 |         0 |   0 |   0 |         0 |        0 |     0 |
| Bismuth                                     |    0 |        38 |         0 |   0 |   0 |         0 |        0 |     0 |
| Cadmium                                     |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Calcium                                     |    0 |        38 |         0 |   0 |   0 |         0 |        0 |     0 |
| Chromium (total)                            |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Cobalt                                      |    0 |       113 |         0 |   0 |   0 |         0 |        0 |     0 |
| Copper                                      |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Iron                                        |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Lead                                        |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Lithium                                     |    0 |        38 |         0 |   0 |   0 |         0 |        0 |     0 |
| Magnesium                                   |    0 |        38 |         0 |   0 |   0 |         0 |        0 |     0 |
| Manganese                                   |    0 |       165 |         0 |   0 |   0 |         0 |        0 |     0 |
| Mercury                                     |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Molybdenum                                  |    0 |       113 |         0 |   0 |   0 |         0 |        0 |     0 |
| Nickel                                      |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Potassium                                   |    0 |        38 |         0 |   0 |   0 |         0 |        0 |     0 |
| Ruthenium                                   |    0 |        38 |         0 |   0 |   0 |         0 |        0 |     0 |
| Selenium                                    |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Silver                                      |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Sodium                                      |    0 |        38 |         0 |   0 |   0 |         0 |        0 |     0 |
| Strontium                                   |    0 |       113 |         0 |   0 |   0 |         0 |        0 |     0 |
| Tellurium                                   |    0 |        38 |         0 |   0 |   0 |         0 |        0 |     0 |
| Thallium                                    |    0 |       113 |         0 |   0 |   0 |         0 |        0 |     0 |
| Tin                                         |    0 |       165 |         0 |   0 |   0 |         0 |        0 |     0 |
| Titanium                                    |    0 |        75 |         0 |   0 |   0 |         0 |        0 |     0 |
| Vanadium                                    |    0 |       113 |         0 |   0 |   0 |         0 |        0 |     0 |
| Zinc                                        |    0 |       230 |         0 |   0 |   0 |         0 |        0 |     0 |
| Butyltin (mono+di+tri)                      |    0 |         0 |       123 |   0 |   0 |         0 |        0 |     0 |
| Dibutyltin                                  |    0 |         0 |       123 |   0 |   0 |         0 |        0 |     0 |
| Monobutyltin                                |    0 |         0 |       123 |   0 |   0 |         0 |        0 |     0 |
| Tetrabutyltin                               |    0 |         0 |        59 |   0 |   0 |         0 |        0 |     0 |
| Tributyltin                                 |    0 |         0 |       123 |   0 |   0 |         0 |        0 |     0 |
| 1-Methylnaphthalene                         |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| 1-Methylphenanthrene                        |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| 1,1-Biphenyl                                |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| 2-Methylnaphthalene                         |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| 2,3,5-Trimethylnaphthalene                  |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| 2,6-Dimethylnaphthalene                     |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Acenaphthene                                |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Acenaphthylene                              |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Anthracene                                  |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Benzo(a)anthracene                          |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Benzo(a)pyrene                              |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Benzo(b)fluoranthene                        |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Benzo(e)pyrene                              |    0 |         0 |         0 | 147 |   0 |         0 |        0 |     0 |
| Benzo(g,h,i)perylene                        |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Benzo(k)fluoranthene                        |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Chrysene                                    |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Dibenz(a,h)anthracene                       |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Dibenzothiophene                            |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Fluoranthene                                |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Fluorene                                    |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Indeno(1,2,3-cd)pyrene                      |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Naphthalene                                 |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| PAHs (High MW 13)                           |    0 |         0 |         0 | 143 |   0 |         0 |        0 |     0 |
| PAHs (High MW)                              |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| PAHs (Low MW 9)                             |    0 |         0 |         0 | 143 |   0 |         0 |        0 |     0 |
| PAHs (Low MW)                               |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| PAHs (total 22)                             |    0 |         0 |         0 | 143 |   0 |         0 |        0 |     0 |
| PAHs (total)                                |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Perylene                                    |    0 |         0 |         0 | 147 |   0 |         0 |        0 |     0 |
| Phenanthrene                                |    0 |         0 |         0 | 147 |   0 |         0 |        0 |     0 |
| Pyrene                                      |    0 |         0 |         0 | 225 |   0 |         0 |        0 |     0 |
| Aroclor-1254                                |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| Aroclor-1260                                |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-003 (4-CB)                              |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-007 (2,4-DiCB)                          |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-008 (2,4’-DiCB)                         |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-015 (4,4’-DiCB)                         |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-016+032                                 |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-018 (2,2’,5-TrCB)                       |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-022 (2,3,4’-TrCB)                       |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-024 (2,3,6-TrCB)                        |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-025 (2,3’,4-TrCB)                       |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-026 (2,3’,5-TrCB)                       |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-028 (2,4,4’-TrCB)                       |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-028+031                                 |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-029 (2,4,5-TrCB)                        |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-031 (2,4’,5-TrCB)                       |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-033 (2,3’,4’-TrCB)                      |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-037 (3,4,4’-TrCB)                       |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-037+042                                 |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-040 (2,2’,3,3’-TeCB)                    |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-041+064                                 |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-044 (2,2’,3,5’-TeCB)                    |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-045 (2,2’,3,6-TeCB)                     |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-046 (2,2’,3,6’-TeCB)                    |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-047+048                                 |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-049 (2,2’,4,5’-TeCB)                    |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-050 (2,2’,4,6-TeCB)                     |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-052 (2,2’,5,5’-TeCB)                    |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-056+060                                 |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-066 (2,3’,4,4’-TeCB)                    |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-070 (2,3’,4’,5-TeCB)                    |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-074 (2,4,4’,5-TeCB)                     |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-077 (3,3’,4,4’-TeCB)                    |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-081 (3,4,4’,5-TeCB)                     |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-082 (2,2’,3,3’,4-PeCB)                  |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-083 (2,2’,3,3’,5-PeCB)                  |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-084 (2,2’,3,3’,6-PeCB)                  |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-085 (2,2’,3,4,4’-PeCB)                  |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-087 (2,2’,3,4,5’-PeCB)                  |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-088 (2,2’,3,4,6-PeCB)                   |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-092 (2,2’,3,5,5’-PeCB)                  |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-095 (2,2’,3,5’,6-PeCB)                  |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-097 (2,2’,3,4’,5’-PeCB)                 |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-099 (2,2’,4,4’,5-PeCB)                  |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-100 (2,2’,4,4’,6-PeCB)                  |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-101 (2,2’,4,5,5’-PeCB)                  |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-105 (2,3,3’,4,4’-PeCB)                  |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-110 (2,3,3’,4’,6-PeCB)                  |    0 |         0 |         0 |   0 |  82 |         0 |        0 |     0 |
| PCB-114 (2,3,4,4’,5-PeCB)                   |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-118 (2,3’,4,4’,5-PeCB)                  |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-119 (2,3’,4,4’,6-PeCB)                  |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-123 (2,3’,4,4’,5’-PeCB)                 |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-126 (3,3’,4,4’,5-PeCB)                  |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-128 (2,2’,3,3’,4,4’-HxCB)               |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-129 (2,2’,3,3’,4,5-HxCB)                |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-132+168                                 |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-136 (2,2’,3,3’,6,6’-HxCB)               |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-137 (2,2’,3,4,4’,5-HxCB)                |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-138 (2,2’,3,4,4’,5’-HxCB)               |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-141 (2,2’,3,4,5,5’-HxCB)                |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-146 (2,2’,3,4’,5,5’-HxCB)               |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-149 (2,2’,3,4’,5’,6-HxCB)               |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-151 (2,2’,3,5,5’,6-HxCB)                |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-153 (2,2’,4,4’,5,5’-HxCB)               |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-156 (2,3,3’,4,4’,5-HxCB)                |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-157 (2,3,3’,4,4’,5’-HxCB)               |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-158 (2,3,3’,4,4’,6-HxCB)                |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-167 (2,3’,4,4’,5,5’-HxCB)               |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-169 (3,3’,4,4’,5,5’-HxCB)               |    0 |         0 |         0 |   0 |  75 |         0 |        0 |     0 |
| PCB-170 (2,2’,3,3’,4,4’,5-HpCB)             |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-172 (2,2’,3,3’,4,5,5’-HpCB)             |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-174 (2,2’,3,3’,4,5,6’-HpCB)             |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-177 (2,2’,3,3’,4,5’,6’-HpCB)            |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-178 (2,2’,3,3’,5,5’,6-HpCB)             |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-180 (2,2’,3,4,4’,5,5’-HpCB)             |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-183 (2,2’,3,4,4’,5’,6-HpCB)             |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-185 (2,2’,3,4,5,5’,6-HpCB)              |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-187 (2,2’,3,4’,5,5’,6-HpCB)             |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-188 (2,2’,3,4’,5,6,6’-HpCB)             |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-189 (2,3,3’,4,4’,5,5’-HpCB)             |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-191 (2,3,3’,4,4’,5’,6-HpCB)             |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-194 (2,2’,3,3’,4,4’,5,5’-OcCB)          |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-195 (2,2’,3,3’,4,4’,5,6-OcCB)           |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-196 (2,2’,3,3’,4,4’,5,6’-OcCB)          |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-199+200                                 |    0 |         0 |         0 |   0 |  33 |         0 |        0 |     0 |
| PCB-200 (2,2’,3,3’,4,5,6,6’-OcCB)           |    0 |         0 |         0 |   0 | 107 |         0 |        0 |     0 |
| PCB-201 (2,2’,3,3’,4,5’,6,6’-OcCB)          |    0 |         0 |         0 |   0 | 140 |         0 |        0 |     0 |
| PCB-205 (2,3,3’,4,4’,5,5’,6-OcCB)           |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCB-206 (2,2’,3,3’,4,4’,5,5’,6-NoCB)        |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCB-209 (DeCB)                              |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCBCL6 (unidentified PCBs)                  |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCBs (total)                                |    0 |         0 |         0 |   0 | 223 |         0 |        0 |     0 |
| PCBTRI2 (unidentified PCBs)                 |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| PCBTRI4 (unidentified PCBs)                 |    0 |         0 |         0 |   0 |  65 |         0 |        0 |     0 |
| 1,1-Dichloro-2,2-bis(4-ethylphenyl) ethane  |    0 |         0 |         0 |   0 |   0 |        75 |        0 |     0 |
| 2,4’-DDD                                    |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| 2,4’-DDE                                    |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| 2,4’-DDT                                    |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| 4,4’-DDD                                    |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| 4,4’-DDE                                    |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| 4,4’-DDT                                    |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Aldrin                                      |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| alpha-BHC                                   |    0 |         0 |         0 |   0 |   0 |       140 |        0 |     0 |
| alpha-Chlordane                             |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| beta-BHC                                    |    0 |         0 |         0 |   0 |   0 |       140 |        0 |     0 |
| Chlordane (total)                           |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| cis-Nonachlor                               |    0 |         0 |         0 |   0 |   0 |       140 |        0 |     0 |
| DDT+DDE+DDD (sum)                           |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| delta-BHC                                   |    0 |         0 |         0 |   0 |   0 |       140 |        0 |     0 |
| Dieldrin                                    |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Endosulfan I                                |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Endosulfan II                               |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Endosulfan sulfate                          |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Endrin                                      |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Endrin aldehyde                             |    0 |         0 |         0 |   0 |   0 |       140 |        0 |     0 |
| Endrin ketone                               |    0 |         0 |         0 |   0 |   0 |        75 |        0 |     0 |
| gamma-BHC                                   |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| gamma-Chlordane                             |    0 |         0 |         0 |   0 |   0 |       140 |        0 |     0 |
| Heptachlor                                  |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Heptachlor epoxide                          |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Hexachlorobenzene                           |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Methoxychlor                                |    0 |         0 |         0 |   0 |   0 |        75 |        0 |     0 |
| Mirex                                       |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Oxychlordane                                |    0 |         0 |         0 |   0 |   0 |       140 |        0 |     0 |
| Pesticides (total)                          |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Toxaphene                                   |    0 |         0 |         0 |   0 |   0 |       141 |        0 |     0 |
| trans-Nonachlor                             |    0 |         0 |         0 |   0 |   0 |       223 |        0 |     0 |
| Organic Carbon (total)                      |    0 |         0 |         0 |   0 |   0 |         0 |      285 |     0 |
| Percent Clay                                |    0 |         0 |         0 |   0 |   0 |         0 |      183 |     0 |
| Percent Coarse Sand                         |    0 |         0 |         0 |   0 |   0 |         0 |       75 |     0 |
| Percent Fine Sand                           |    0 |         0 |         0 |   0 |   0 |         0 |       75 |     0 |
| Percent Gravel                              |    0 |         0 |         0 |   0 |   0 |         0 |      118 |     0 |
| Percent Medium Sand                         |    0 |         0 |         0 |   0 |   0 |         0 |       75 |     0 |
| Percent Pebbles and Shells                  |    0 |         0 |         0 |   0 |   0 |         0 |       75 |     0 |
| Percent Sand                                |    0 |         0 |         0 |   0 |   0 |         0 |      214 |     0 |
| Percent Silt                                |    0 |         0 |         0 |   0 |   0 |         0 |      183 |     0 |
| Percent Silt and Clay                       |    0 |         0 |         0 |   0 |   0 |         0 |      214 |     0 |
| Percent Total Grain Size                    |    0 |         0 |         0 |   0 |   0 |         0 |       75 |     0 |
| Percent Very Coarse Sand                    |    0 |         0 |         0 |   0 |   0 |         0 |       72 |     0 |
| Percent Very Fine Sand                      |    0 |         0 |         0 |   0 |   0 |         0 |       72 |     0 |
| Solids (total)                              |    0 |         0 |         0 |   0 |   0 |         0 |       75 |     0 |
| Uranium-238                                 |    0 |         0 |         0 |   0 |   0 |         0 |        0 |    38 |

``` r
rm(tmp)
```

Note that parameters include totals and sums. But, confusingly, not all
sums included here were (re)calculated in the Access database Sums
table. Some must pre-date Ramboll’s work.

# Units

``` r
xtabs(~Era + Parameter_Group  + Units, data = sedtox_data)
```

    ## , , Units = %
    ## 
    ##        Parameter_Group
    ## Era     CDDF Inorganic Organotin  PAH  PCB Pesticide Physical Radio
    ##   1990s    0         0         0    0    0         0      323     0
    ##   2000s    0         0         0    0    0         0      330     0
    ##   2010s    0         0         0    0    0         0     1063     0
    ## 
    ## , , Units = µg/g dry
    ## 
    ##        Parameter_Group
    ## Era     CDDF Inorganic Organotin  PAH  PCB Pesticide Physical Radio
    ##   1990s    0       715         0    0    0         0        0     0
    ##   2000s    0      1815         0    0    0         0        0    38
    ##   2010s    0      1830         0    0    0         0       75     0
    ## 
    ## , , Units = ng/g dry
    ## 
    ##        Parameter_Group
    ## Era     CDDF Inorganic Organotin  PAH  PCB Pesticide Physical Radio
    ##   1990s    0         0       128 2015 5200      1950        0     0
    ##   2000s    0         0       128 2184 1596      1748        0     0
    ##   2010s    0         0       295 2296 4204      2554        0     0
    ## 
    ## , , Units = ng/kg dry
    ## 
    ##        Parameter_Group
    ## Era     CDDF Inorganic Organotin  PAH  PCB Pesticide Physical Radio
    ##   1990s  690         0         0    0    0         0        0     0
    ##   2000s  736         0         0    0    0         0        0     0
    ##   2010s  459         0         0    0    0         0        0     0

Note that metals are in \(\mu g/g\) or ppm, organics are in \(n g/g\) or
ppb, physical observations in percent.

# Subdivide Data By Parameter Group

It is worth separating by Parameter Group into smaller, more manageable
data sets. We could have done this more compactly with lapply, but by
repeating code, we could specify file names simply and directly.

``` r
diox_data   <- sedtox_data %>% filter(Parameter_Group == 'CDDF')
metals_data <- sedtox_data %>% filter(Parameter_Group == 'Inorganic')
tbt_data    <- sedtox_data %>% filter(Parameter_Group == 'Organotin')
pah_data    <- sedtox_data %>% filter(Parameter_Group == 'PAH')
pcb_data    <- sedtox_data %>% filter(Parameter_Group == 'PCB')
pest_data   <- sedtox_data %>% filter(Parameter_Group == 'Pesticide')
phys_data   <- sedtox_data %>% filter(Parameter_Group == 'Physical')
```

The challenge we have with many of these smaller data sets will be
extracting meaning in a simple way. We have far too many parameters to
include in the Chapter in anything but a confusing way. This particular
data presentation includes both individual parameters and sums, so we
need to extract information one parameter at a time, not do blind sums
by parameter groups.

Unfortunately, the list of parameters actually obscures which are sums
and which are field measurements of mixtures.

# Check for Other Sums and Totals in the Data

## The “Official” Totals

Ramboll Calculated the following “Official” Totals, which

1.  Should be in the dataset, and  
2.  Are easy to remove because we have a complete list.

<!-- end list -->

``` r
xtabs(~ Parameter + Parameter_Group , data = sums_data)
```

    ##                         Parameter_Group
    ## Parameter                CDDF Organotin PAH PCB Pesticide
    ##   Butyltin (mono+di+tri)    0       123   0   0         0
    ##   CDD/CDF (total)          79         0   0   0         0
    ##   Chlordane (total)         0         0   0   0       223
    ##   DDT+DDE+DDD (sum)         0         0   0   0       223
    ##   PAHs (High MW)            0         0 225   0         0
    ##   PAHs (Low MW)             0         0 225   0         0
    ##   PAHs (total)              0         0 225   0         0
    ##   PCBs (total)              0         0   0 223         0
    ##   Pesticides (total)        0         0   0   0       223

## Checking All Data Subsets

We look at each of the data subsets to identify “totals” and sums that
are in the original data, but do not represent original observations. We
do this simply by examining the list of parameters in each group.

### Butyl Tin Totals

``` r
xtabs(~ Parameter + Parameter_Group , data = tbt_data)
```

    ##                         Parameter_Group
    ## Parameter                Organotin
    ##   Butyltin (mono+di+tri)       123
    ##   Dibutyltin                   123
    ##   Monobutyltin                 123
    ##   Tetrabutyltin                 59
    ##   Tributyltin                  123

Only a single total, which is the “official” total recalculated by
Ramboll.

### Dioxin Totals

``` r
xtabs(~ Parameter + Parameter_Group , data = diox_data)
```

    ##                                              Parameter_Group
    ## Parameter                                     CDDF
    ##   1,2,3,4,6,7,8-Heptachlorodibenzo-p-dioxin     79
    ##   1,2,3,4,6,7,8-Heptachlorodibenzofuran         79
    ##   1,2,3,4,7,8-Hexachlorodibenzo-p-dioxin        79
    ##   1,2,3,4,7,8-Hexachlorodibenzofuran            79
    ##   1,2,3,4,7,8,9-Heptachlorodibenzofuran         79
    ##   1,2,3,6,7,8-Hexachlorodibenzo-p-dioxin        79
    ##   1,2,3,6,7,8-Hexachlorodibenzofuran            79
    ##   1,2,3,7,8-Pentachlorodibenzo-p-dioxin         79
    ##   1,2,3,7,8-Pentachlorodibenzofuran             79
    ##   1,2,3,7,8,9-Hexachlorodibenzo-p-dioxin        79
    ##   1,2,3,7,8,9-Hexachlorodibenzofuran            79
    ##   2,3,4,6,7,8-Hexachlorodibenzofuran            79
    ##   2,3,4,7,8-Pentachlorodibenzofuran             79
    ##   2,3,7,8-Tetrachlorodibenzo-p-dioxin           79
    ##   2,3,7,8-Tetrachlorodibenzofuran               79
    ##   CDD/CDF (total)                               79
    ##   HpCDD (total)                                 17
    ##   HpCDF (total)                                 17
    ##   HxCDD (total)                                 17
    ##   HxCDF (total)                                 17
    ##   Octachlorodibenzo-p-dioxin                    79
    ##   Octachlorodibenzofuran                        79
    ##   PCB-077 (3,3',4,4'-TeCB) (dioxin-like)        62
    ##   PCB-126 (3,3',4,4',5-PeCB) (dioxin-like)      62
    ##   PCB-169 (3,3',4,4',5,5'-HxCB) (dioxin-like)   62
    ##   PeCDD (total)                                 17
    ##   PeCDF (total)                                 17
    ##   TCDD (total)                                  17
    ##   TCDF (total)                                  79
    ##   TEQ CDD/CDF                                   79

The only “Official” Dioxin-related total is \* CDD/CDF (total)

A small number of samples (17) include other values flagged as (total).
\* HpCDD (total) \* HpCDF (total) \* HxCDD (total) \* HxCDF (total) \*
PeCDD (total) \* PeCDF (total) \* TCDD (total) \* TCDF (total)

Three other items turn up in 62 samples. They are not labeled as totals,
and they refer to specific tetra- penta- and hexa- chlorobiphenyls.  
(note that these are listed under `Parameter_Group == "CDDF"`, NOT
`Parameter_Group == "PCB"`.) \* PCB-077 (3,3’,4,4’-TeCB) (dioxin-like)  
\* PCB-126 (3,3’,4,4’,5-PeCB) (dioxin-like)  
\* PCB-169 (3,3’,4,4’,5,5’-HxCB) (dioxin-like)

It is probably not an accident that the sum of 62 + 17 = 79, which
appears to be the total number of samples.

Reviewing the full list of parameters, None of these appear to be sums.

HpCDD, for example, refers to Heptachlorodibenzo-para-dioxin, which
occurs in a number of isomers, which are chemically difficult to
distinguish in samples.  
The “(total)” here appears to refer to a total of multiple (analytically
indistinguishable) isomers. NO other HpCDD appears in the parameter
list, so this can not be a sum of other anylates.

A quick check suggest we are looking at differences in reporting from
different years. The 17 samples with `HpCDD (total)` and its cousins are
all from 2010.

``` r
selected <- c('HpCDD (total)', "PCB-077 (3,3',4,4'-TeCB) (dioxin-like)")
diox_data %>%
  select(Sample_ID, Sample_Year, Method, Parameter, Result) %>%
  filter(Parameter %in% selected) %>%
  pivot_wider(c(Sample_ID, Sample_Year, Method), names_from = Parameter, values_from = Result)
```

    ## # A tibble: 79 x 5
    ##    Sample_ID Sample_Year Method `PCB-077 (3,3',4,4'-TeCB) (diox~ `HpCDD (total)`
    ##    <chr>           <dbl> <chr>                             <dbl>           <dbl>
    ##  1 1994.CS04        1994 <NA>                               9.61              NA
    ##  2 1994.EB04        1994 <NA>                              21                 NA
    ##  3 1994.EB06        1994 <NA>                              72.4               NA
    ##  4 1994.EB07        1994 <NA>                              44.2               NA
    ##  5 1994.EB09        1994 <NA>                              72.6               NA
    ##  6 1994.IB01        1994 <NA>                             108                 NA
    ##  7 1994.IB02        1994 <NA>                              98.4               NA
    ##  8 1994.IB03        1994 <NA>                              80.3               NA
    ##  9 1994.IB04        1994 <NA>                              51.9               NA
    ## 10 1994.IB06        1994 <NA>                              41.2               NA
    ## # ... with 69 more rows

### Conclusion

Remaining “Totals” are NOT sums that we should remove.

## Metals Data Totals

``` r
xtabs(~ Parameter + Parameter_Group , data = metals_data)
```

    ##                   Parameter_Group
    ## Parameter          Inorganic
    ##   Aluminum               165
    ##   Antimony               165
    ##   Arsenic                230
    ##   Barium                 113
    ##   Beryllium              113
    ##   Bismuth                 38
    ##   Cadmium                230
    ##   Calcium                 38
    ##   Chromium (total)       230
    ##   Cobalt                 113
    ##   Copper                 230
    ##   Iron                   230
    ##   Lead                   230
    ##   Lithium                 38
    ##   Magnesium               38
    ##   Manganese              165
    ##   Mercury                230
    ##   Molybdenum             113
    ##   Nickel                 230
    ##   Potassium               38
    ##   Ruthenium               38
    ##   Selenium               230
    ##   Silver                 230
    ##   Sodium                  38
    ##   Strontium              113
    ##   Tellurium               38
    ##   Thallium               113
    ##   Tin                    165
    ##   Titanium                75
    ##   Vanadium               113
    ##   Zinc                   230

No Totals. "Chromium (total) refers to measuring all oxidation states of
chromium, not a sum of other observations.

## PAH Data Totals

``` r
xtabs(~ Parameter + Parameter_Group , data = pah_data)
```

    ##                             Parameter_Group
    ## Parameter                    PAH
    ##   1-Methylnaphthalene        225
    ##   1-Methylphenanthrene       225
    ##   1,1-Biphenyl               225
    ##   2-Methylnaphthalene        225
    ##   2,3,5-Trimethylnaphthalene 225
    ##   2,6-Dimethylnaphthalene    225
    ##   Acenaphthene               225
    ##   Acenaphthylene             225
    ##   Anthracene                 225
    ##   Benzo(a)anthracene         225
    ##   Benzo(a)pyrene             225
    ##   Benzo(b)fluoranthene       225
    ##   Benzo(e)pyrene             147
    ##   Benzo(g,h,i)perylene       225
    ##   Benzo(k)fluoranthene       225
    ##   Chrysene                   225
    ##   Dibenz(a,h)anthracene      225
    ##   Dibenzothiophene           225
    ##   Fluoranthene               225
    ##   Fluorene                   225
    ##   Indeno(1,2,3-cd)pyrene     225
    ##   Naphthalene                225
    ##   PAHs (High MW 13)          143
    ##   PAHs (High MW)             225
    ##   PAHs (Low MW 9)            143
    ##   PAHs (Low MW)              225
    ##   PAHs (total 22)            143
    ##   PAHs (total)               225
    ##   Perylene                   147
    ##   Phenanthrene               147
    ##   Pyrene                     225

So under Pahs, we have six totals,

The following PAH Totals are NOT included in the “official” totals, and
should be removed from the data:

  - PAHs (High MW 13)  
  - PAHs (Low MW 9)  
  - PAHs (total 22)

## PCB Data Totals

``` r
xtabs(~ Parameter + Parameter_Group , data = pcb_data)
```

    ##                                       Parameter_Group
    ## Parameter                              PCB
    ##   Aroclor-1254                          65
    ##   Aroclor-1260                          65
    ##   PCB-003 (4-CB)                        75
    ##   PCB-007 (2,4-DiCB)                    65
    ##   PCB-008 (2,4'-DiCB)                  223
    ##   PCB-015 (4,4'-DiCB)                   65
    ##   PCB-016+032                           65
    ##   PCB-018 (2,2',5-TrCB)                223
    ##   PCB-022 (2,3,4'-TrCB)                 65
    ##   PCB-024 (2,3,6-TrCB)                  65
    ##   PCB-025 (2,3',4-TrCB)                 65
    ##   PCB-026 (2,3',5-TrCB)                 65
    ##   PCB-028 (2,4,4'-TrCB)                223
    ##   PCB-028+031                           65
    ##   PCB-029 (2,4,5-TrCB)                  65
    ##   PCB-031 (2,4',5-TrCB)                 75
    ##   PCB-033 (2,3',4'-TrCB)               140
    ##   PCB-037 (3,4,4'-TrCB)                 75
    ##   PCB-037+042                           65
    ##   PCB-040 (2,2',3,3'-TeCB)              65
    ##   PCB-041+064                           65
    ##   PCB-044 (2,2',3,5'-TeCB)             223
    ##   PCB-045 (2,2',3,6-TeCB)               65
    ##   PCB-046 (2,2',3,6'-TeCB)              65
    ##   PCB-047+048                           65
    ##   PCB-049 (2,2',4,5'-TeCB)             140
    ##   PCB-050 (2,2',4,6-TeCB)               65
    ##   PCB-052 (2,2',5,5'-TeCB)             223
    ##   PCB-056+060                          140
    ##   PCB-066 (2,3',4,4'-TeCB)             223
    ##   PCB-070 (2,3',4',5-TeCB)             140
    ##   PCB-074 (2,4,4',5-TeCB)              140
    ##   PCB-077 (3,3',4,4'-TeCB)             223
    ##   PCB-081 (3,4,4',5-TeCB)               75
    ##   PCB-082 (2,2',3,3',4-PeCB)            65
    ##   PCB-083 (2,2',3,3',5-PeCB)            65
    ##   PCB-084 (2,2',3,3',6-PeCB)            65
    ##   PCB-085 (2,2',3,4,4'-PeCB)            65
    ##   PCB-087 (2,2',3,4,5'-PeCB)           140
    ##   PCB-088 (2,2',3,4,6-PeCB)             65
    ##   PCB-092 (2,2',3,5,5'-PeCB)            65
    ##   PCB-095 (2,2',3,5',6-PeCB)            75
    ##   PCB-097 (2,2',3,4',5'-PeCB)          140
    ##   PCB-099 (2,2',4,4',5-PeCB)           140
    ##   PCB-100 (2,2',4,4',6-PeCB)            65
    ##   PCB-101 (2,2',4,5,5'-PeCB)           223
    ##   PCB-105 (2,3,3',4,4'-PeCB)           223
    ##   PCB-110 (2,3,3',4',6-PeCB)            82
    ##   PCB-114 (2,3,4,4',5-PeCB)             75
    ##   PCB-118 (2,3',4,4',5-PeCB)           223
    ##   PCB-119 (2,3',4,4',6-PeCB)            75
    ##   PCB-123 (2,3',4,4',5'-PeCB)           75
    ##   PCB-126 (3,3',4,4',5-PeCB)           223
    ##   PCB-128 (2,2',3,3',4,4'-HxCB)        223
    ##   PCB-129 (2,2',3,3',4,5-HxCB)          65
    ##   PCB-132+168                           75
    ##   PCB-136 (2,2',3,3',6,6'-HxCB)         65
    ##   PCB-137 (2,2',3,4,4',5-HxCB)          65
    ##   PCB-138 (2,2',3,4,4',5'-HxCB)        223
    ##   PCB-141 (2,2',3,4,5,5'-HxCB)         140
    ##   PCB-146 (2,2',3,4',5,5'-HxCB)         65
    ##   PCB-149 (2,2',3,4',5',6-HxCB)        140
    ##   PCB-151 (2,2',3,5,5',6-HxCB)         140
    ##   PCB-153 (2,2',4,4',5,5'-HxCB)        223
    ##   PCB-156 (2,3,3',4,4',5-HxCB)          75
    ##   PCB-157 (2,3,3',4,4',5'-HxCB)         75
    ##   PCB-158 (2,3,3',4,4',6-HxCB)         140
    ##   PCB-167 (2,3',4,4',5,5'-HxCB)        140
    ##   PCB-169 (3,3',4,4',5,5'-HxCB)         75
    ##   PCB-170 (2,2',3,3',4,4',5-HpCB)      223
    ##   PCB-172 (2,2',3,3',4,5,5'-HpCB)       65
    ##   PCB-174 (2,2',3,3',4,5,6'-HpCB)      140
    ##   PCB-177 (2,2',3,3',4,5',6'-HpCB)     140
    ##   PCB-178 (2,2',3,3',5,5',6-HpCB)       65
    ##   PCB-180 (2,2',3,4,4',5,5'-HpCB)      223
    ##   PCB-183 (2,2',3,4,4',5',6-HpCB)      140
    ##   PCB-185 (2,2',3,4,5,5',6-HpCB)        65
    ##   PCB-187 (2,2',3,4',5,5',6-HpCB)      223
    ##   PCB-188 (2,2',3,4',5,6,6'-HpCB)       65
    ##   PCB-189 (2,3,3',4,4',5,5'-HpCB)      140
    ##   PCB-191 (2,3,3',4,4',5',6-HpCB)       65
    ##   PCB-194 (2,2',3,3',4,4',5,5'-OcCB)   140
    ##   PCB-195 (2,2',3,3',4,4',5,6-OcCB)    223
    ##   PCB-196 (2,2',3,3',4,4',5,6'-OcCB)    65
    ##   PCB-199+200                           33
    ##   PCB-200 (2,2',3,3',4,5,6,6'-OcCB)    107
    ##   PCB-201 (2,2',3,3',4,5',6,6'-OcCB)   140
    ##   PCB-205 (2,3,3',4,4',5,5',6-OcCB)     65
    ##   PCB-206 (2,2',3,3',4,4',5,5',6-NoCB) 223
    ##   PCB-209 (DeCB)                       223
    ##   PCBCL6 (unidentified PCBs)            65
    ##   PCBs (total)                         223
    ##   PCBTRI2 (unidentified PCBs)           65
    ##   PCBTRI4 (unidentified PCBs)           65

Here the only total is the “official” total, so we are probably O.K. It
is not clear how we should handle the three “(unidentified PCBs)”.

## Pesticides Data Totals

``` r
xtabs(~ Parameter + Parameter_Group , data = pest_data)
```

    ##                                             Parameter_Group
    ## Parameter                                    Pesticide
    ##   1,1-Dichloro-2,2-bis(4-ethylphenyl) ethane        75
    ##   2,4'-DDD                                         223
    ##   2,4'-DDE                                         223
    ##   2,4'-DDT                                         223
    ##   4,4'-DDD                                         223
    ##   4,4'-DDE                                         223
    ##   4,4'-DDT                                         223
    ##   Aldrin                                           223
    ##   alpha-BHC                                        140
    ##   alpha-Chlordane                                  223
    ##   beta-BHC                                         140
    ##   Chlordane (total)                                223
    ##   cis-Nonachlor                                    140
    ##   DDT+DDE+DDD (sum)                                223
    ##   delta-BHC                                        140
    ##   Dieldrin                                         223
    ##   Endosulfan I                                     223
    ##   Endosulfan II                                    223
    ##   Endosulfan sulfate                               223
    ##   Endrin                                           223
    ##   Endrin aldehyde                                  140
    ##   Endrin ketone                                     75
    ##   gamma-BHC                                        223
    ##   gamma-Chlordane                                  140
    ##   Heptachlor                                       223
    ##   Heptachlor epoxide                               223
    ##   Hexachlorobenzene                                223
    ##   Methoxychlor                                      75
    ##   Mirex                                            223
    ##   Oxychlordane                                     140
    ##   Pesticides (total)                               223
    ##   Toxaphene                                        141
    ##   trans-Nonachlor                                  223

The three totals here  
\* Chlordane (total) \* DDT+DDE+DDD (sum)  
\* Pesticides (total) are all “official” totals.

## Physical Data Totals

``` r
xtabs(~ Parameter + Sample_Year , data = phys_data)
```

    ##                             Sample_Year
    ## Parameter                    1991 2000 2001 2002 2010 2011
    ##   Organic Carbon (total)       63   27   20   18   73   84
    ##   Percent Clay                 65   20    6   17   33   42
    ##   Percent Coarse Sand           0    0    0    0   33   42
    ##   Percent Fine Sand             0    0    0    0   33   42
    ##   Percent Gravel                0   20    7   18   31   42
    ##   Percent Medium Sand           0    0    0    0   33   42
    ##   Percent Pebbles and Shells    0    0    0    0   33   42
    ##   Percent Sand                 65   27   22   18   40   42
    ##   Percent Silt                 65   20    6   17   33   42
    ##   Percent Silt and Clay        65   27   22   18   40   42
    ##   Percent Total Grain Size      0    0    0    0   33   42
    ##   Percent Very Coarse Sand      0    0    0    0   30   42
    ##   Percent Very Fine Sand        0    0    0    0   30   42
    ##   Solids (total)                0    0    0    0   33   42

None of the parameters here with ‘Total’ in the names are sums, however,
some of the OTHER values represent sums calculated before the data was
submitted to Ramboll.

  - **Percent Silt and Clay** is (almost always) a sum of separate silt
    and clay values. However, seven 2000 samples, one 2002 sample and
    seven 2010 samples only have “Percent Sand” and “Percent Silt and
    Clay” values.
  - Sand was subdivided into more precise sand fractions only in 2010
    and 2011. **Percent Sand** in those two years was USUALLY a SUM of
    the other sand fractions , except for seven samples from the NCCA
    from 2010.
  - **Percent Total Grain Size** was a QA/QC check testing whether the
    sum of all grain size data sums to one (or close to it).

Note that “Percent Gravel” was inconsistently reported, making it
unclear whether all other grain size fractions are comparable, at least
for the handful of sites with high gravel.

These sums probably need to be removed and recalculated, without
removing the data associated with samples where we have no finer-scale
data.

### Checks on sums and totals

#### Silt And Clay Totals

Silt Clay is consistently EITHER the sum of other silt and clay values
or a measuremetn where silt and clay values are not also available.

``` r
phys_data %>%
  filter(Parameter != 'Organic Carbon (total)') %>%
  filter(Parameter != 'Solids (total)') %>%
  group_by(Sample_ID) %>%
  summarize(year = first(Sample_Year),
         sand = sum(grepl('Sand', Parameter) * Result),
         silt = sum((grepl('Silt', Parameter) & ! grepl('Clay', Parameter)) * 
                      Result),
         clay = sum((grepl('Clay', Parameter) & ! grepl('Silt', Parameter)) *
                      Result),
         siltclay = sum((Parameter == "Percent Silt and Clay") *
                          Result),
         test = siltclay - (silt + clay) < 0.01 | (silt == 0 & clay == 0),
         .groups = 'drop')
```

    ## # A tibble: 213 x 7
    ##    Sample_ID  year  sand  silt  clay siltclay test 
    ##    <chr>     <dbl> <dbl> <dbl> <dbl>    <dbl> <lgl>
    ##  1 1991.CS01  1991  88.1  10.7   1.2     11.9 TRUE 
    ##  2 1991.CS02  1991  87.1  10.6   2.4     13   TRUE 
    ##  3 1991.CS03  1991  84.1  13     2.9     15.9 TRUE 
    ##  4 1991.CS04  1991  29.9  32.6  37.5     70.1 TRUE 
    ##  5 1991.CS05  1991  82.4  12.7   4.9     17.6 TRUE 
    ##  6 1991.CS06  1991  65.9  21.6  12.5     34.1 TRUE 
    ##  7 1991.CS07  1991  89.5  10.2   0.3     10.5 TRUE 
    ##  8 1991.EB01  1991  34.7  26    39.4     65.4 TRUE 
    ##  9 1991.EB02  1991  33.2  29.4  37.4     66.8 TRUE 
    ## 10 1991.EB03  1991  25.1  42.1  32.8     74.9 TRUE 
    ## # ... with 203 more rows

We can selectively remove “Silt + Clay” values that are sums of other
variables as follows:

``` r
phys_data %>%
  filter(Parameter != 'Organic Carbon (total)') %>%
  filter(Parameter != 'Solids (total)') %>%
  group_by(Sample_ID) %>%
  mutate(silt = sum((grepl('Silt', Parameter) & ! grepl('Clay', Parameter)) * 
                      Result),
         clay = sum((grepl('Clay', Parameter) & ! grepl('Silt', Parameter)) *
                      Result),
         test = ! (silt == 0 & clay == 0)) %>%
  ungroup() %>%
  filter( ! (Parameter == 'Percent Silt and Clay' & test)) %>%
  select(-silt, -clay, -test)
```

    ## # A tibble: 1,248 x 26
    ##    Project Region Location Substation Lat_N Long_W Sample_ID Sample_Year
    ##    <chr>   <chr>  <chr>    <chr>      <dbl>  <dbl> <chr>           <dbl>
    ##  1 1991-2~ Cape ~ L.CS01   S.CS01      43.7  -69.9 1991.CS01        1991
    ##  2 1991-2~ Cape ~ L.CS01   S.CS01      43.7  -69.9 1991.CS01        1991
    ##  3 1991-2~ Cape ~ L.CS01   S.CS01      43.7  -69.9 1991.CS01        1991
    ##  4 1991-2~ Cape ~ L.CS02   S.CS02      43.7  -69.9 1991.CS02        1991
    ##  5 1991-2~ Cape ~ L.CS02   S.CS02      43.7  -69.9 1991.CS02        1991
    ##  6 1991-2~ Cape ~ L.CS02   S.CS02      43.7  -69.9 1991.CS02        1991
    ##  7 1991-2~ Cape ~ L.CS03   S.CS03      43.7  -69.9 1991.CS03        1991
    ##  8 1991-2~ Cape ~ L.CS03   S.CS03      43.7  -69.9 1991.CS03        1991
    ##  9 1991-2~ Cape ~ L.CS03   S.CS03      43.7  -69.9 1991.CS03        1991
    ## 10 1991-2~ Cape ~ L.CS04Q  S.CS04      43.7  -69.9 1991.CS04        1991
    ## # ... with 1,238 more rows, and 18 more variables: Sample_Date <dttm>,
    ## #   Replicate <dbl>, Sample_Type <chr>, Matrix <chr>, Parameter_Group <chr>,
    ## #   Parameter <chr>, CASRN <chr>, Result <dbl>, MDL <dbl>, RL <dbl>,
    ## #   Units <chr>, Detect <chr>, Det_Flag <dbl>, Qualifier <chr>, `QA
    ## #   Qualifier` <chr>, Reportable_Result <chr>, Method <chr>, Era <fct>

#### Sand Totals

Here’s a quick way to find the samples from 2010 and 2011 we do NOT want
to change.

``` r
phys_data %>%
  filter(Sample_Year > 2009, grepl("NCCA", Sample_ID)) %>%
  arrange(Parameter)
```

    ## # A tibble: 21 x 26
    ##    Project Region Location Substation Lat_N Long_W Sample_ID Sample_Year
    ##    <chr>   <chr>  <chr>    <chr>      <dbl>  <dbl> <chr>           <dbl>
    ##  1 NCCA D~ Outer~ L.1015   <NA>        43.7  -70.1 NCCA10-1~        2010
    ##  2 NCCA D~ East ~ L.1013   <NA>        43.9  -69.9 NCCA10-1~        2010
    ##  3 NCCA D~ East ~ L.1016   <NA>        43.7  -70.0 NCCA10-1~        2010
    ##  4 NCCA D~ West ~ L.1019   <NA>        43.8  -70.0 NCCA10-1~        2010
    ##  5 NCCA D~ Outer~ L.1017   <NA>        43.6  -70.2 NCCA10-1~        2010
    ##  6 NCCA D~ Inner~ L.1021   <NA>        43.7  -70.2 NCCA10-1~        2010
    ##  7 NCCA D~ Outer~ L.1022   <NA>        43.7  -70.1 NCCA10-1~        2010
    ##  8 NCCA D~ Outer~ L.1015   <NA>        43.7  -70.1 NCCA10-1~        2010
    ##  9 NCCA D~ East ~ L.1013   <NA>        43.9  -69.9 NCCA10-1~        2010
    ## 10 NCCA D~ East ~ L.1016   <NA>        43.7  -70.0 NCCA10-1~        2010
    ## # ... with 11 more rows, and 18 more variables: Sample_Date <dttm>,
    ## #   Replicate <dbl>, Sample_Type <chr>, Matrix <chr>, Parameter_Group <chr>,
    ## #   Parameter <chr>, CASRN <chr>, Result <dbl>, MDL <dbl>, RL <dbl>,
    ## #   Units <chr>, Detect <chr>, Det_Flag <dbl>, Qualifier <chr>, `QA
    ## #   Qualifier` <chr>, Reportable_Result <chr>, Method <chr>, Era <fct>

So we can invert that and add a filter for the parameter, and find the
data rows we want to delete.

``` r
phys_data %>%
  filter(Sample_Year > 2009,
         ! grepl("NCCA", Sample_ID),
         Parameter == 'Percent Sand')
```

    ## # A tibble: 75 x 26
    ##    Project Region Location Substation Lat_N Long_W Sample_ID Sample_Year
    ##    <chr>   <chr>  <chr>    <chr>      <dbl>  <dbl> <chr>           <dbl>
    ##  1 Casco ~ Outer~ L.OB05   S.OB05      43.7  -70.1 CBEP2010~        2010
    ##  2 Casco ~ Cape ~ L.CS03T  S.CS03T     43.7  -69.9 CBEP2010~        2011
    ##  3 Casco ~ Inner~ L.IB03T  S.IB03      43.7  -70.2 CBEP2010~        2010
    ##  4 Casco ~ West ~ L.WB03   S.WB03      43.8  -70.0 CBEP2010~        2010
    ##  5 Casco ~ Inner~ L.SW01T  S.SW01      43.6  -70.3 CBEP2010~        2010
    ##  6 Casco ~ Outer~ L.OB05   S.OB05      43.7  -70.1 CBEP2010~        2010
    ##  7 Casco ~ West ~ L.WB04T  S.WB04      43.8  -70.0 CBEP2010~        2011
    ##  8 Casco ~ Inner~ L.SW02   S.SW02      43.6  -70.3 CBEP2010~        2010
    ##  9 Casco ~ East ~ L.EB04   S.EB04      43.7  -69.9 CBEP2010~        2011
    ## 10 Casco ~ Inner~ L.IB02   S.IB02      43.7  -70.2 CBEP2010~        2010
    ## # ... with 65 more rows, and 18 more variables: Sample_Date <dttm>,
    ## #   Replicate <dbl>, Sample_Type <chr>, Matrix <chr>, Parameter_Group <chr>,
    ## #   Parameter <chr>, CASRN <chr>, Result <dbl>, MDL <dbl>, RL <dbl>,
    ## #   Units <chr>, Detect <chr>, Det_Flag <dbl>, Qualifier <chr>, `QA
    ## #   Qualifier` <chr>, Reportable_Result <chr>, Method <chr>, Era <fct>

Note that we found 75 values, exactly what we would expect based on
sample frequencies in the cross-tabs, above.

Lets check those sites to make sure the “Percent Sand” value really is a
sum.

``` r
phys_data %>%
  filter(Sample_Year > 2009) %>%
  group_by(Sample_ID) %>%
  summarize(sand = sum((grepl('Sand', Parameter) & ! Parameter == 'Percent Sand') * 
                      Result),
         total = sum((Parameter == 'Percent Sand') * 
                      Result),
         test =  abs(total - sand) < 0.001,
         .groups = 'drop')
```

    ## # A tibble: 82 x 4
    ##    Sample_ID      sand total test 
    ##    <chr>         <dbl> <dbl> <lgl>
    ##  1 CBEP2010-CS01  99.2  99.2 TRUE 
    ##  2 CBEP2010-CS02  17    17   TRUE 
    ##  3 CBEP2010-CS03  97.3  97.3 TRUE 
    ##  4 CBEP2010-CS04  50.6  50.6 TRUE 
    ##  5 CBEP2010-CS05  95.6  95.6 TRUE 
    ##  6 CBEP2010-CS06  52.4  52.4 TRUE 
    ##  7 CBEP2010-CS07  89.3  89.3 TRUE 
    ##  8 CBEP2010-EB01  18.8  18.8 TRUE 
    ##  9 CBEP2010-EB02  34.6  34.6 TRUE 
    ## 10 CBEP2010-EB03  21.4  21.4 TRUE 
    ## # ... with 72 more rows

It is.

## Total Grain Size

Percent Total Grain Size was not calculated for all samples, and appears
to Include the Gravel Fraction. It should be removed from the data as
inconsistent and possibly misleading.

``` r
unique(phys_data$Parameter)
```

    ##  [1] "Percent Sand"               "Percent Silt"              
    ##  [3] "Percent Clay"               "Percent Silt and Clay"     
    ##  [5] "Percent Gravel"             "Organic Carbon (total)"    
    ##  [7] "Solids (total)"             "Percent Coarse Sand"       
    ##  [9] "Percent Fine Sand"          "Percent Medium Sand"       
    ## [11] "Percent Pebbles and Shells" "Percent Very Coarse Sand"  
    ## [13] "Percent Very Fine Sand"     "Percent Total Grain Size"

``` r
phys_data %>%
  filter(Parameter != 'Organic Carbon (total)') %>%
  filter(Parameter != 'Solids (total)') %>%
  group_by(Sample_ID) %>%
  summarize(year = first(Sample_Year),
         sand = sum(grepl('Sand', Parameter) * Result),
         siltclay = sum((Parameter == "Percent Silt and Clay") *
                          Result),
         gravel = sum((Parameter == 'Percent Gravel') * Result),
         total = sum((Parameter == 'Percent Total Grain Size') * Result),
         test = abs(total + gravel - 100) < 0.1,
         .groups = 'drop')
```

    ## # A tibble: 213 x 7
    ##    Sample_ID  year  sand siltclay gravel total test 
    ##    <chr>     <dbl> <dbl>    <dbl>  <dbl> <dbl> <lgl>
    ##  1 1991.CS01  1991  88.1     11.9      0     0 FALSE
    ##  2 1991.CS02  1991  87.1     13        0     0 FALSE
    ##  3 1991.CS03  1991  84.1     15.9      0     0 FALSE
    ##  4 1991.CS04  1991  29.9     70.1      0     0 FALSE
    ##  5 1991.CS05  1991  82.4     17.6      0     0 FALSE
    ##  6 1991.CS06  1991  65.9     34.1      0     0 FALSE
    ##  7 1991.CS07  1991  89.5     10.5      0     0 FALSE
    ##  8 1991.EB01  1991  34.7     65.4      0     0 FALSE
    ##  9 1991.EB02  1991  33.2     66.8      0     0 FALSE
    ## 10 1991.EB03  1991  25.1     74.9      0     0 FALSE
    ## # ... with 203 more rows

``` r
phys_data %>%
  filter(Parameter != 'Organic Carbon (total)') %>%
  filter(Parameter != 'Solids (total)') %>%
  filter(Sample_ID == 'CBEP2010-WB08')
```

    ## # A tibble: 12 x 26
    ##    Project Region Location Substation Lat_N Long_W Sample_ID Sample_Year
    ##    <chr>   <chr>  <chr>    <chr>      <dbl>  <dbl> <chr>           <dbl>
    ##  1 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ##  2 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ##  3 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ##  4 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ##  5 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ##  6 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ##  7 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ##  8 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ##  9 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ## 10 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ## 11 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ## 12 Casco ~ West ~ L.WB08T  S.WB08      43.7  -70.0 CBEP2010~        2011
    ## # ... with 18 more variables: Sample_Date <dttm>, Replicate <dbl>,
    ## #   Sample_Type <chr>, Matrix <chr>, Parameter_Group <chr>, Parameter <chr>,
    ## #   CASRN <chr>, Result <dbl>, MDL <dbl>, RL <dbl>, Units <chr>, Detect <chr>,
    ## #   Det_Flag <dbl>, Qualifier <chr>, `QA Qualifier` <chr>,
    ## #   Reportable_Result <chr>, Method <chr>, Era <fct>

The inconsistent reporting of the Gravel Fraction makes it impossible to
be sure how to interpret the sand and silt/sand fractions for all
samples. Best practice here may be to recalculate those values based on
the gravel-free fraction. We test that idea here.

``` r
phys_data %>%
  group_by(Sample_ID, Replicate) %>%
  summarize(year = first(Sample_Year),
         total = sum(grepl('Percent Total Grain Size', Parameter) * Result),
         total2 = sum(Result) - 
           sum((Parameter == 'Percent Gravel') * Result ) -
           sum((Parameter == 'Percent Silt') * Result ) -
           sum((Parameter == 'Percent Clay') * Result ),
         gravel = sum((Parameter == 'Percent Gravel') * Result),
         .groups = 'drop')
```

    ## # A tibble: 215 x 6
    ##    Sample_ID Replicate  year total total2 gravel
    ##    <chr>         <dbl> <dbl> <dbl>  <dbl>  <dbl>
    ##  1 1991.CS01         0  1991     0   100.      0
    ##  2 1991.CS02         0  1991     0   100.      0
    ##  3 1991.CS03         0  1991     0   100.      0
    ##  4 1991.CS04         0  1991     0   103.      0
    ##  5 1991.CS05         0  1991     0   100.      0
    ##  6 1991.CS06         0  1991     0   101.      0
    ##  7 1991.CS07         0  1991     0   100.      0
    ##  8 1991.EB01         0  1991     0   102.      0
    ##  9 1991.EB02         0  1991     0   102.      0
    ## 10 1991.EB03         0  1991     0   102.      0
    ## # ... with 205 more rows

That points out that ZERO values here are “really” NAs.

# Examine Structure of Sums Data

``` r
xtabs(~Era + Parameter_Group, data = sums_data)
```

    ##        Parameter_Group
    ## Era     CDDF Organotin PAH PCB Pesticide
    ##   1990s   30        32 195  65       195
    ##   2000s   32        32 234  76       228
    ##   2010s   17        59 246  82       246

# Parameters and Parameter Groups

``` r
the_order <- as.numeric(factor(sums_data$Parameter_Group))

tmp <- sums_data %>%
  mutate(Parameter = factor(Parameter,
                            levels = levels(fct_reorder(Parameter, the_order))))
                                                 
knitr::kable(xtabs(~Parameter + Parameter_Group, data = tmp))
```

|                        | CDDF | Organotin | PAH | PCB | Pesticide |
| :--------------------- | ---: | --------: | --: | --: | --------: |
| CDD/CDF (total)        |   79 |         0 |   0 |   0 |         0 |
| Butyltin (mono+di+tri) |    0 |       123 |   0 |   0 |         0 |
| PAHs (High MW)         |    0 |         0 | 225 |   0 |         0 |
| PAHs (Low MW)          |    0 |         0 | 225 |   0 |         0 |
| PAHs (total)           |    0 |         0 | 225 |   0 |         0 |
| PCBs (total)           |    0 |         0 |   0 | 223 |         0 |
| Chlordane (total)      |    0 |         0 |   0 |   0 |       223 |
| DDT+DDE+DDD (sum)      |    0 |         0 |   0 |   0 |       223 |
| Pesticides (total)     |    0 |         0 |   0 |   0 |       223 |

``` r
rm(tmp)
```

Those look pretty convenient for comparison purposes.

**The metadata is not entirely clear, but it appears (from examining the
Access database) that these sums omit non-detects, effectively equating
non-detects to zero. That is inconsistent with how we handled
non-detects in several other toxics data sets, where we have been using
maximum likelihood estimators of expected values.**

# Units

``` r
xtabs(~Era + Parameter_Group  + Units, data = sums_data)
```

    ## , , Units = ng/g dry
    ## 
    ##        Parameter_Group
    ## Era     CDDF Organotin PAH PCB Pesticide
    ##   1990s    0        32 195  65       195
    ##   2000s    0        32 234  76       228
    ##   2010s    0        59 246  82       246
    ## 
    ## , , Units = ng/kg dry
    ## 
    ##        Parameter_Group
    ## Era     CDDF Organotin PAH PCB Pesticide
    ##   1990s   30         0   0   0         0
    ##   2000s   32         0   0   0         0
    ##   2010s   17         0   0   0         0

Ramboll Standardized units in the Access database, so, MOST sums are
expressed in ng/g dry weight (\~ ppb).

The Dioxins and Furans are expressed in ng/kg, or pg/g or approximately
parts per trilion. There are no Squirts for Dioxins and Furans. Instead,
Ramboll ALSO expressed them in TEQ – Tox equivalents. Toxic equivalents
provide a way to estimate the cumulative toxic effect of a mixture of
related chemicals by eighting each compound by its relative toxic
effect, compared to some reference compound (conventionally TCDD).

# Load SQuiRTs Data

We compare those totals to the SQuiRTs, when availalbe. Ramboll
Environ’s Access database again contains a table derived from the
SQuiRTs, containing ERL (Effects Range Low) and ERM (Effects Range
Medium) values for selected parameters.

``` r
fn <- "SQuiRTs.xlsx"
squirts_data <- read_excel(file.path(sibling,fn), 
    col_types = c("skip", "text", "text", 
        "text", "numeric", "numeric", "text", 
        "numeric", "numeric", "text", "skip"))
```

## Units

We need to ensure that we are using similar units for analysis and
comparison to the SQuiRTs. The SQuiRTs Table contains two different sets
of values, with different units. We want to check both.

``` r
xtabs(~ Param_Norm + Units_org , data = squirts_data)
```

    ##                        Units_org
    ## Param_Norm              ppb
    ##   2-Methylnaphthalene     1
    ##   4,4'-DDD                1
    ##   4,4'-DDE                1
    ##   4,4'-DDT                1
    ##   Acenaphthene            1
    ##   Acenaphthylene          1
    ##   Anthracene              1
    ##   Arsenic                 1
    ##   Benzo(a)anthracene      1
    ##   Benzo(a)pyrene          1
    ##   Cadmium                 1
    ##   Chlordane (total)       1
    ##   Chromium (total)        1
    ##   Chrysene                1
    ##   Copper                  1
    ##   DDT+DDE+DDD (sum)       1
    ##   Dibenz(a,h)anthracene   1
    ##   Dieldrin                1
    ##   Fluoranthene            1
    ##   Fluorene                1
    ##   Lead                    1
    ##   Mercury                 1
    ##   Naphthalene             1
    ##   Nickel                  1
    ##   PAHs (High MW)          1
    ##   PAHs (Low MW)           1
    ##   PAHs (total)            1
    ##   PCBs (total)            1
    ##   Phenanthrene            1
    ##   Pyrene                  1
    ##   Silver                  1
    ##   Zinc                    1

``` r
xtabs(~ Param_Norm + Units_fin, data = squirts_data)
```

    ##                        Units_fin
    ## Param_Norm              µg/g dry ng/g dry
    ##   2-Methylnaphthalene          0        1
    ##   4,4'-DDD                     0        1
    ##   4,4'-DDE                     0        1
    ##   4,4'-DDT                     0        1
    ##   Acenaphthene                 0        1
    ##   Acenaphthylene               0        1
    ##   Anthracene                   0        1
    ##   Arsenic                      1        0
    ##   Benzo(a)anthracene           0        1
    ##   Benzo(a)pyrene               0        1
    ##   Cadmium                      1        0
    ##   Chlordane (total)            0        1
    ##   Chromium (total)             1        0
    ##   Chrysene                     0        1
    ##   Copper                       1        0
    ##   DDT+DDE+DDD (sum)            0        1
    ##   Dibenz(a,h)anthracene        0        1
    ##   Dieldrin                     0        1
    ##   Fluoranthene                 0        1
    ##   Fluorene                     0        1
    ##   Lead                         1        0
    ##   Mercury                      1        0
    ##   Naphthalene                  0        1
    ##   Nickel                       1        0
    ##   PAHs (High MW)               0        1
    ##   PAHs (Low MW)                0        1
    ##   PAHs (total)                 0        1
    ##   PCBs (total)                 0        1
    ##   Phenanthrene                 0        1
    ##   Pyrene                       0        1
    ##   Silver                       1        0
    ##   Zinc                         1        0

So, SQuiRTs are originally in ppb, but Ramboll Environ converted units
to µg/g dry or ng/g dry, consistent with the way concentrations were
also reported.
