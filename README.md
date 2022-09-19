# S2_random_forest
Land cover classification using Random Forest Classifier. In order to use this script you need to have shapefile of the boundary and sample of landcover (with point geometry) uploaded to Google Earth Engine.

# How to use this script :
Please kindly set up and configure some of these variables :

## 1. Date Range
Date range selelection for Sentinel 2 aquisition date <br>
`(on line 49-50)` <br>
`var START_DATE = ee.Date('2022-04-01');` <br>
`var END_DATE = ee.Date('2022-09-01');`

## 2. Center Coordinates for the Map
Please plug in the coordinates values from console to the script. This need to be done in order for the map window to show the location of the AO <br>
`(on line 233)` <br>
`linker.get(0).setCenter(111.93846248507332, -7.597410151981024, 12);`

Cheers
