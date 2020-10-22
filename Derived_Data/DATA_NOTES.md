#  Shapefiles
Shapefiles were produced as follows:

1.  Data from the `Locations.xls` file was processed in the R Markdown file 
    "Prepare_GIS_data.Rmd", producing two CSV files:  
    *   'locations.csv' and  
    *   'substations.csv'  

2.  Those files were imported into ArcGIS, as tables.

3.  We generated event layers from those data tables (based on UTM coordinates,
    Zone 19 N, NAD 1983), and

4.  Eexported the resulting layers as shapefiles.
