# Data Source
All data is from the official data transmission for Ramboll Environ to CBEP in
2017, at the completion of Ramboll Environ's analysis of Casco Bay Sediment
chemistry data.

All data was derived from an Access database prepared by Ramboll Environ.
The full database, called “Casco Bay Sediment Data 1991_2001.accdb” and dated
February 1, 2017, is available from CBEP upon request.

##  Sediment Toxics Data
The file **"Casco Bay Sediment Data Query.xlsx"** is the output of the
"Results_Env_withSums" query in that database, which is described in the
metadata as:

>  This is the final data query that contains all of the standardized sediment
   data from the Casco Bay sediment dataset (1991-2011). 

## Data on Sums of Certain Toxic Compounds
The File **"Sums_Query.xlsx"** is the output of the "Sums" query, which is
described in the metadata as:

>  This query calculates the chemical sums (i.e., PAHs) in the Casco Bay
   sediment data. It feeds the final Res_Env_withSums query.

We use it as a convenient source for the sums also included in the previous
file. To simplify public presentation of complex results, *State of Casco Bay*
will emphasize various sums and totals, rather than report on dozens of toxic
constituents. It is often simpler to begin with the smaller data set.

**The metadata is not entirely clear, but it appears (from examining the Access
database) that these sums omit non-detects, effectively equating non-detects to
zero. That is inconsistent with how we  handled non-detects in several other
toxics data sets, where we have been using maximum likelihood estimators of
expected values.**

## Documentation of Components of Those Sums
The file **"SumGroups.xlsx"** is an excel copy of the internal "SumGroups" table 
from the Access database. The metadata includes the following:

>  This table identifies which individual chemicals are included in the chemical
   sums calculated in the database.
   
We include it here for its value as metadata, explaining what compounds are
included in the various sums and totals.  

Because of differences in methodologies used over the years, many sums and
totals are based on relatively small subsets of constituent chemicals, and thus
may underestimate results that would have come from more comprehensive surveys.

### Number of Chemical Constituents in "Totals"  
Analytic Total    | Number of Compounds
------------------|------------------------
CDD/CDF (total)*  |  17
PAHs (total)      |  13
PAHs (High MW)    |   6
PAHs (Low MW)     |   7
PCBs (Total)      |  88
Pesticides (total)|  20
Chlordane (total) |   3
DDT+DDE+DDD (sum) |   6

* A combination of dioxins and furans.

## Data on Sampling Locations
The file **"Locations.xlsx"** is an excel copy of the internal "locations' table
from the Ramboll Access database. Locations are in decimal degrees and UTM
coordinates, based on WGS 1983. The table received minimal description in the 
associated metadata, saying only: 

>  This table lists all of the sediment sample locations in Casco Bay from
   1991-2011.

Location data was originally developed by CBEP, so we can add the following
information:

Places where samples were collected are presented as:

1.  Stations -- nominal sampling locations consistent used for sample planning 
    across eras. For a variety of reasons, distances between samples associated
    with a single Station can sometimes be substantial.  
    
2.  Locations -- Nominal Locations associated with stations.  Locations
    aggregate multiple sampling locations into approximate common locations. All
    samples assigned to the  same Location were collected within a radius of
    approximately 100 m.

3.  Substations -- These provide more detailed estimates of location, with all
    observations associated with a single substation collected within a distance 
    of approximately 10 meters of one another.

## Information on Reference Levels of COntaminants in Sediments
The file **"SQuiRTs.xlsx"** contains ERL and ERM values taken from the Ramboll
Environ Access database. In particular it contains ERL (Effects Range Low) and
ERM (Effects Range Medium) values for selected parameters.  These values were,
in turn, derived from NOAA's SQuiRT tables.  Additional details on the SQuiRT
tables are available in the archives containing data on analysis of toxic
compounds in Porland Harbor sediments, available at
https://github.com/ccb60/PortlandHarborToxics .
