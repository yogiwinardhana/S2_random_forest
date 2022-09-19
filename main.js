// ######################################################################################################
// ######################################################################################################
//                                    ### SENTINEL 2 FOR DUE DILIGENCE  ###
// ######################################################################################################
/* 
   Metadata:
    Title: Sentinel 2 For Due Diligence
    Base: -
    Update: - Add SRTM Data
            - Add buffer boundary
            - Add single S2SR imageryl
    Project: Due Diligence
    Author: I Made Y.W.
*/

var centroid = function(feature) {
  return feature.centroid()
}

var centroid_table = table.map(centroid)

//-- Get latitude and longitude of a centroid
var table_latlng = centroid_table.geometry().coordinates()
var table_lat = table_latlng.get(0)
var table_lng = table_latlng.get(1)
print("Latitude : ", table_lat)
print("Longitude : ", table_lng)

//-- Do buffer 1 km for polygon
var table_buffer = function(feature) {
  return feature.buffer(1000);
};

var table_buffer = table.map(table_buffer);

// ######################################################################################################
//                                    ### Add Sentinel 2 Imagery  ###
// ######################################################################################################
//-- Add image collection for copernicus s2 surface reflectance
var s2Sr = ee.ImageCollection('COPERNICUS/S2_SR')

//-- Add image collection for copernicus s2 surface cloud probability
var s2Clouds = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY');

//-- Add start date to sort image
//-- tahun n
//-- 2021-05-01 // 2021-11-15

var START_DATE = ee.Date('2022-04-01');
var END_DATE = ee.Date('2022-09-01');

//-- Add single S2SR
var s2SrSingle = ee.ImageCollection('COPERNICUS/S2_SR')
                  .filterBounds(table)
                  .filterDate(START_DATE, END_DATE)
                  .sort('CLOUDY_PIXEL_PERCENTAGE')
                  .first()
                  .clip(table_buffer)
                  
print("S2 Surface Reflectance Single", s2SrSingle)

//-- Set threshold for MAX_CLOUD_PROBABILITY
var MAX_CLOUD_PROBABILITY = 10;

//-- Set region using imported shapefile
var region = table

//-- Function to maskcloud
function maskClouds(img) {
  var clouds = ee.Image(img.get('cloud_mask')).select('probability');
  var isNotCloud = clouds.lt(MAX_CLOUD_PROBABILITY);
  return img.updateMask(isNotCloud);
}

//-- The masks for the 10m bands sometimes do not exclude bad data at
//-- scene edges, so we apply masks from the 20m and 60m bands as well.
//-- Example asset that needs this operation:
//-- COPERNICUS/S2_CLOUD_PROBABILITY/20190301T000239_20190301T000238_T55GDP
function maskEdges(s2_img) {
  return s2_img.updateMask(
      s2_img.select('B8A').mask().updateMask(s2_img.select('B9').mask()));
}

//-- Filter input collections by desired data range and region.
var criteria = ee.Filter.and(ee.Filter.bounds(region), 
                             ee.Filter.date(START_DATE, END_DATE));
s2Sr = s2Sr.filter(criteria).map(maskEdges);
print("S2 Surface Reflectance", s2Sr)

s2Clouds = s2Clouds.filter(criteria);
print("S2 Cloud Probability", s2Clouds)

//-- Join S2_SR with cloud probability dataset to add cloud mask.
var s2SrWithCloudMask = ee.Join.saveFirst('cloud_mask').apply({
  primary: s2Sr,
  secondary: s2Clouds,
  condition:
      ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'})
});

//-- Get median for s2SrWithCloudMask
var s2CloudMasked = ee.ImageCollection(s2SrWithCloudMask).map(maskClouds).median();

//-- Clip S2CloudMasked imagery
var s2CloudMaskedClip = s2CloudMasked.clip(table_buffer)

//---------------------------------8 Vegetation indices
// https://www.l3harrisgeospatial.com/docs/broadbandgreenness.html#Differen

//1. NDVI = (NIR -Red) / (NIR + Red)
var ndvi = s2CloudMaskedClip.expression(
  ' (NIR - RED) / (NIR + RED) ', {
    'NIR' : s2CloudMaskedClip.select('B8'),
    'RED' : s2CloudMaskedClip.select('B4')
  }).rename('NDVI');

//2.  GRVI = (NIR /  Green)
var grvi = s2CloudMaskedClip.expression(
  '  NIR / GREEN ', {
    'NIR' : s2CloudMaskedClip.select('B8'),
    'GREEN' : s2CloudMaskedClip.select('B3')
  }).rename('GRVI');

//3. OSAVI = (NIR - Red) / (NIR + Red + 0.16)
var osavi = s2CloudMaskedClip.expression(
  ' (NIR - RED) / (NIR + RED + 0.16) ', {
    'NIR' : s2CloudMaskedClip.select('B8'),
    'RED' : s2CloudMaskedClip.select('B4')
  }).rename('OSAVI');

//4. EVI = 2.5 * ( (NIR - Red) / (NIR + 6 * Red - 7.5 * Blue + 1) )
var evi = s2CloudMaskedClip.expression(
  ' 2.5 * ( (NIR - RED) / (NIR + 6 * oRED - 7.5 * BLUE + 1) ) ', {
    'NIR' : s2CloudMaskedClip.select('B8'),
    'RED' : s2CloudMaskedClip.select('B4'),
    'BLUE' : s2CloudMaskedClip.select('B2')
  }).rename('EVI')

//5. LAI = (3.618 * EVI - 0.118)
var lai =  evi.expression(
  ' (3.618 * EVI - 0.118) ', {
    'EVI' : evi.select('EVI')
  }).rename('LAI')

//6. GNDVI = (NIR - Green) / (NIR + Green)
var gndvi = s2CloudMaskedClip.expression(
  ' (NIR - GREEN) / (NIR + GREEN) ', {
    'NIR' : s2CloudMaskedClip.select('B8'),
    'GREEN' : s2CloudMaskedClip.select('B3')
  }).rename('GNDVI')

//7. RDVI = ( (NIR - Red) / (NIR + Red)^ 0.5 ) 
var rdvi = s2CloudMaskedClip.expression(
  ' (NIR - RED) / (NIR + RED ) ^ 0.5 ', {
    'NIR' : s2CloudMaskedClip.select('B8'),
    'RED' : s2CloudMaskedClip.select('B4')
  }).rename('RDVI')

//8. DVI = ( NIR - Red)
var dvi = s2CloudMaskedClip.expression(
  ' NIR - RED', {
    'NIR' : s2CloudMaskedClip.select('B8'),
    'RED' : s2CloudMaskedClip.select('B4')
  }).rename('DVI')
  
//9. RBI = ( Red / Blue)
// var rbi = s2CloudMaskedClip.expression(
//   ' (RED / BLUE)', {
//   ' RED' : s2CloudMaskedClip.select('B4'),
//   ' BLUE': s2CloudMaskedClip.select('B2')
// }).rename('RBI')


// NDBI cause much inaccuracy
// //10. NDBI
var ndbi = s2CloudMaskedClip.expression(
      '((SWIR - NIR) / (SWIR + NIR))', {
        'SWIR' : s2CloudMaskedClip.select('B12'),
        'NIR'  : s2CloudMaskedClip.select('B8')
      }
).rename('NDBI')

// ######################################################################################################
//                                    ### Set Sentinel 2 Visualization ###
// ######################################################################################################
var rgbVis = {min: 0, max: 4095, bands: ['B11', 'B8', 'B2']};
var rgbVis_natural = {min: 0, max: 2000, bands: ['B4', 'B3', 'B2'], gamma: 1};


// ######################################################################################################
//                                    ### LINKING UI ###
// ######################################################################################################
//-- Linking map
var linkedMap = ui.Map();

//-- Visualization parameter
var rgbVis = {min: 0, max: 4095, bands: ['B11', 'B8', 'B2']};
var rgbVis_redege = {min: 0, max: 4095, bands: ['B11', 'B8', 'B2']};
var rgbVis_natural = {min: 0, max: 2000, bands: ['B4', 'B3', 'B2'], gamma: 1};

//-- Add first layer which is our clipped raster
Map.addLayer(s2CloudMaskedClip, rgbVis, "S2 SR B11,8,2 Clip Masked at " + MAX_CLOUD_PROBABILITY + '%', false);
//Map.addLayer(s2CloudMaskedClip, rgbVis, "S2 SR B11,8,2 Clip Masked at " + MAX_CLOUD_PROBABILITY + '%', true);

// Add second layer and link it to our first image
linkedMap.addLayer(s2CloudMaskedClip, rgbVis_natural, "S2 Natural Color")
linkedMap.addLayer(s2SrSingle, rgbVis_natural, "S2 Natural Color Single")

// Style for feature
var shown = true;
var not_shown = false 
var opacity = 0.8; 
var nameLayer = 'map'; 
var visParamsFeature = {color: 'red'}; 

linkedMap.addLayer(table_buffer, visParamsFeature, 'Boundary Buffer 1 km', not_shown, opacity);
linkedMap.addLayer(table, visParamsFeature, 'Boundary Non Buffer', not_shown, opacity );

// Set map linker
var linker = ui.Map.Linker([ui.root.widgets().get(0), linkedMap]);

// Set split panel 
var splitPanel = ui.SplitPanel({
  firstPanel: linker.get(0),
  secondPanel: linker.get(1),
  orientation: 'horizontal',
  wipe: true,
  style: {stretch: 'both'}
});
ui.root.widgets().reset([splitPanel]);

// Set center for linker
linker.get(0).setCenter(111.93846248507332, -7.597410151981024, 12);
//linker.get(0).setCenter(table_lng.getInfo(), table_lat.getInfo(), 12);

//-----------Supervised
//-----------Supervised
// Select the bands for training
var bands = ['NDVI',
            'GNDVI',
            'GRVI',
            'DVI',
            'OSAVI',
            //'NDBI',
            //'EVI',
            //'LAI',
            //'RDVI'
            //'RBI'];
            'B2',
            //'B3',,
            //'B4',
            //'B5',
            //'B6',
            //'B7',
            'B8',
            //'B9',
            'B11'];
            //'B12'];

// var S2_bands_classification = s2CloudMaskedClip.select(['B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B11', 'B12' ])
var S2_bands_classification = s2CloudMaskedClip.select(['B2', 'B8', 'B11'])


var veg_composite = ndvi.addBands([ndbi, 
                                  gndvi, 
                                  grvi,
                                  dvi, 
                                  osavi
                                  //ndbi
                                  //evi,
                                  //lai,
                                  //rdvi,
                                  //rbi
                                  ])

var veg_composite = veg_composite.addBands(S2_bands_classification)

print(veg_composite)

// Add sample
// var sample = FOREST
// .merge(AGRI_YOUNG)
// .merge(AGRI_OLD)
// .merge(OPEN_AGRI)
// .merge(PLANTATION)
// .merge(OLD_SCRUB)
// .merge(SETTLEMENT)
// .merge(OPEN_MINING)
// .merge(WATER_BODY)

var sample = table2

// // Generate ramdom column
var gcp = sample.randomColumn()

print("GCP", gcp)

// // Get training and validation
var trainingGCP = gcp.filter(ee.Filter.lt('random', 0.7));
var validationGCP = gcp.filter(ee.Filter.gte('random', 0.7));

print("Number of GCPs", sample.size())
print("Number of training GCP", trainingGCP.size())
print("Number of validation GCP", validationGCP.size())

// //Sample the input imagery to get a FeatureCollection of training data.
var training = veg_composite.sampleRegions({
  collection: trainingGCP, 
  properties: ['LC_CODE'], 
  scale: 10,
  tileScale:16
});

print(training)

// //Making a Random Forest classifier and training it. -----------------------------------NEED MORE ADJUSDMENTS !!!
var classifier= ee.Classifier.smileRandomForest(10).train({
  features: training, 
  classProperty: 'LC_CODE',
  inputProperties: bands
});

// // //Classifying the input imagery
var classification= veg_composite.classify(classifier)

print(classification)

// // // // Define a palette for the Land Use classification.
var palette = [
  '0b5c11', // FOREST
  'acd6a2', //  AGRI_YOUNG
  '009999', // AGRI_OLD
  '95583e', // OPEN AGRI
  'a4c22b', //  PLANTATION
  '99ff99', // OLD SCRUB
  'ff1a1a', // SETTLEMENT
  'ffff99', // OPEN MINING
  '4386ff'// WATER BODY
];

linkedMap.addLayer(classification, {min:1, max:9, palette:palette}, "Supervised")

//************************************************************************** 
// Accuracy Assessment
//************************************************************************** 

// // Use classification map to assess accuracy using the validation fraction
// // of the overall training set created above.
var test = classification.sampleRegions({
  collection: validationGCP,
  properties: ['LC_CODE'],
  tileScale: 16,
  scale: 10,
});

var testConfusionMatrix = test.errorMatrix('LC_CODE')
// // Printing of confusion matrix may time out. Alternatively, you can export it as CSV
print('Confusion Matrix', testConfusionMatrix);
print('Test Accuracy', testConfusionMatrix.accuracy());

//******Part 4.1: Post processing classification******
// https://courses.spatialthoughts.com/end-to-end-gee.html#post-processing-classification-results
// count patch sizes
var patchsize = classification.connectedPixelCount(40, false);

// run a majority filter
var filtered = classification.focal_mode({
    radius: 30,
    kernelType: 'square',
    units: 'meters',
}); 
  
// updated image with majority filter where patch size is small
var connectedClassified =  classification.where(patchsize.lt(25),filtered);

////******Part 5:Create a legend******
//////////////////////////////////////

//Set position of panel
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});
 
//Create legend title
var legendTitle = ui.Label({
  value: 'Classification Legend',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});
 
//Add the title to the panel
legend.add(legendTitle);
 
//Create and style 1 row of the legend.
var makeRow = function(color, name) {
 
      var colorBox = ui.Label({
        style: {
          backgroundColor: '#' + color,
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });
      
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};
 
//Identify palette with the legend colors
var palette =['0b5c11', 'acd6a2', '009999', '95583e', 
'a4c22b','99ff99', 'ff1a1a', 'ffff99', '4386ff'];
 

//Identify names within the legend
var names = ['FOREST','AGRI_YOUNG','AGRI_OLD',
            'OPEN_AGRI','PLANTATION','OLD SCRUB', 'SETTLEMENT', 'OPEN MINING', 'WATER_BODY'];
 
//Add color and names
for (var i = 0; i < 9; i++) {
  legend.add(makeRow(palette[i], names[i]));
  }  

//Add legend to map
Map.add(legend);

////******Part 6: Display the Final Land Cover Classification and Provide Export Options******
//////////////////////////////////////////////////////////////////////////////////////////////

//Create palette for the final land cover map classifications
//var palette =['0e721d', '111eff', 'c6d247', '726d17', 'a8ff29','0b090c'];
var urbanPalette = 
'<RasterSymbolizer>' +
' <ColorMap  type="intervals">' +
    '<ColorMapEntry color="#0e721d" quantity="1" label="FOREST"/>' +
    '<ColorMapEntry color="#111eff" quantity="2" label="AGRI_YOUNG"/>' +
    '<ColorMapEntry color="#c6d247" quantity="3" label="AGRI_OLD"/>' +
    '<ColorMapEntry color="#726d17" quantity="4" label="OPEN_AGRI"/>' +
    '<ColorMapEntry color="#3bff11" quantity="5" label="PLANTATION"/>' +
    '<ColorMapEntry color="#0b090c" quantity="6" label="OLD SCRUB"/>'+
    '<ColorMapEntry color="#0b090c" quantity="7" label="SETTLEMENT"/>'+
    '<ColorMapEntry color="#0b090c" quantity="8" label="OPEN_MINING"/>'+
    '<ColorMapEntry color="#0b090c" quantity="9" label="WATER_BODY"/>'+
  '</ColorMap>' +
'</RasterSymbolizer>';

//Mask out impervious surfaces
//var finalmap = classified.blend(masked);

//Add final map to the display
Map.addLayer(classification.sldStyle(urbanPalette), {}, "Land Classification", 0);
Map.addLayer(connectedClassified.sldStyle(urbanPalette), {}, "Land Classification Post Processing", 0);

// -------------DOWNLOAD IMAGE
var download_img = function(image_name, band_arr, description, region, start_date, end_date){
  /*Param
    bands: (array) list of bands
    image_name:  (variable) name of image,
    description: (string) of description,
    region: (variable) AOI that already buffered,
    start_date: (ee.Date()) starting date
    end_date: (ee.Date()) end date
    
  */
  Export.image.toDrive({
    image:image_name.select(band_arr),
    description: START_DATE.format('YYYYMMdd', 'UTC').getInfo() + "_" + END_DATE.format('YYYYMMdd', 'UTC').getInfo() + "_" + description,
    scale:10,
    fileFormat: 'GeoTIFF',
    maxPixels: 10000000000000,
    region: region,
    formatOptions: {
      cloudOptimized: true
    }
  })
}

download_img(ndvi, ['NDVI'], 'NDVI', table_buffer, START_DATE, END_DATE)
download_img(classification, ['classification'], 'random_forest', table_buffer, START_DATE, END_DATE)
download_img(s2CloudMaskedClip, ['B2', 'B3', 'B4', 'B8', 'B11'], 'S2', table_buffer, START_DATE, END_DATE)

// Related links
// https://developers.google.com/earth-engine/apidocs/ui-splitpanel
// https://developers.google.com/earth-engine/apidocs/ui-panel
// https://developers.google.com/earth-engine/apidocs/ui-splitpanel-setsecondpanel
// https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_CLOUD_PROBABILITY#bands
// https://www.linkedin.com/pulse/split-screen-display-satellite-images-gee-khang-vu-
