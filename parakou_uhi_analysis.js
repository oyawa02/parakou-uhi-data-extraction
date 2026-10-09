// ==============================================================================
// GEE COMPLETE SCRIPT - PARAKOU (LST, NDVI, NDBI, ECI) 
// ==============================================================================

// 1. Definition of the Study Area (ROI) via its exact Cloud Asset path
var roi = ee.FeatureCollection('projects/axial-autonomy-464605-s3/assets/Parakou'); 

// Strict and aligned temporal window (Peak dry season)
var startDOY = '-01-20';
var endDOY   = '-02-15';

// 2. Landsat 8 Preprocessing Function (Collection 2, Level 2)
function processLandsat8(image) {
  // Masking clouds and cloud shadows using QA_PIXEL
  var qa = image.select('QA_PIXEL');
  var cloudShadowMask = (1 << 4);
  var cloudsMask = (1 << 3);
  var clearMask = qa.bitwiseAnd(cloudShadowMask).eq(0)
    .and(qa.bitwiseAnd(cloudsMask).eq(0));
  
  // Applying optical and thermal scaling factors
  var optical = image.select(['SR_B4', 'SR_B3', 'SR_B2', 'SR_B5', 'SR_B6'])
    .multiply(0.0000275).add(-0.2)
    .clip(roi);
    
  var thermal = image.select('ST_B10')
    .multiply(0.00341802).add(149.0) // Temperature in Kelvin
    .subtract(273.15)              // Conversion to °C
    .clip(roi)
    .rename('LST');

  // Calculating indices
  var red  = optical.select('SR_B4'); // Red Band
  var nir  = optical.select('SR_B5'); // NIR Band
  var swir = optical.select('SR_B6'); // SWIR1 Band

  var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
  var ndbi = swir.subtract(nir).divide(swir.add(nir)).rename('NDBI');

  return optical.addBands(ndvi).addBands(ndbi).addBands(thermal)
    .updateMask(clearMask)
    .copyProperties(image, ['system:time_start']);
}

// 3. Loading and Composing Collections (2014 vs 2024)
var col2014 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2014' + startDOY, '2014' + endDOY)
  .filter(ee.Filter.lt('CLOUD_COVER', 10))
  .map(processLandsat8);

var col2024 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2024' + startDOY, '2024' + endDOY)
  .filter(ee.Filter.lt('CLOUD_COVER', 10))
  .map(processLandsat8);

print('Landsat Collection 2014:', col2014);
print('Landsat Collection 2024:', col2024);

// Creating median composites
var comp2014 = col2014.median().clip(roi);
var comp2024 = col2024.median().clip(roi);

// Initial raw band extraction
var rawLst2014  = comp2014.select('LST');
var rawLst2024  = comp2024.select('LST');

// INTERANNUAL THERMAL BIAS CORRECTION (Z-Score Normalization over the ROI)
var harmonizeLST = function(img, geometry) {
  var stats = img.reduceRegion({
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    geometry: geometry,
    scale: 30,
    maxPixels: 1e9
  });
  var mean = ee.Number(stats.get('LST_mean'));
  var std  = ee.Number(stats.get('LST_stdDev'));
  return img.subtract(mean).divide(std).multiply(2.5).add(40).rename('LST');
};

var lst2014 = harmonizeLST(rawLst2014, roi);
var lst2024 = harmonizeLST(rawLst2024, roi);

var ndvi2014 = comp2014.select('NDVI').rename('NDVI');
var ndvi2024 = comp2024.select('NDVI').rename('NDVI');

var ndbi2014 = comp2014.select('NDBI').rename('NDBI');
var ndbi2024 = comp2024.select('NDBI').rename('NDBI');

// 4. Robust Percentile Normalization (2% - 98%) for ECI Calculation
var getRobustNorm = function(img, bandName) {
  var singleBand = img.select(bandName);
  var stats = singleBand.reduceRegion({
    reducer: ee.Reducer.percentile([2, 98]),
    geometry: roi,
    scale: 30,
    maxPixels: 1e9
  });
  var p2  = ee.Number(stats.get(bandName + '_p2'));
  var p98 = ee.Number(stats.get(bandName + '_p98'));
  return singleBand.clamp(p2, p98).subtract(p2).divide(p98.subtract(p2));
};

var normNDVI_14 = getRobustNorm(comp2014, 'NDVI');
var normNDVI_24 = getRobustNorm(comp2024, 'NDVI');

var normNDBI_14 = getRobustNorm(comp2014, 'NDBI');
var normNDBI_24 = getRobustNorm(comp2024, 'NDBI');

var normLST_14  = getRobustNorm(lst2014, 'LST');
var normLST_24  = getRobustNorm(lst2024, 'LST');

// Environmental Criticality Index (ECI) Calculation
var eci2014 = normNDBI_14.add(ee.Image(1).subtract(normNDVI_14)).add(normLST_14).divide(3).rename('ECI');
var eci2024 = normNDBI_24.add(ee.Image(1).subtract(normNDVI_24)).add(normLST_24).divide(3).rename('ECI');

// ==============================================================================
// 5. MAP VISUALIZATION AND PARAKOU BOUNDARIES
// ==============================================================================
Map.centerObject(roi, 12);

Map.addLayer(comp2014.select(['SR_B6', 'SR_B5', 'SR_B4']), {min: 0, max: 0.3}, 'Landsat 2014 (False Color)');
Map.addLayer(comp2024.select(['SR_B6', 'SR_B5', 'SR_B4']), {min: 0, max: 0.3}, 'Landsat 2024 (False Color)');
Map.addLayer(eci2014, {min: 0, max: 1, palette: ['blue', 'yellow', 'red']}, 'ECI 2014');
Map.addLayer(eci2024, {min: 0, max: 1, palette: ['blue', 'yellow', 'red']}, 'ECI 2024');

// ADMINISTRATIVE BOUNDARIES OUTLINE
var outline = ee.Image().byte().paint({
  featureCollection: roi,
  color: 1,
  width: 2
});
Map.addLayer(outline, {palette: 'red'}, 'Parakou Boundaries');

// ==============================================================================
// 6. INTERACTIVE CONSOLE CHARTS
// ==============================================================================

// A. Harmonized LST Histogram (°C)
var chartLST = ui.Chart.image.histogram({
  image: ee.Image.cat([lst2014.rename('LST 2014'), lst2024.rename('LST 2024')]),
  region: roi,
  scale: 30,
  maxPixels: 1e9
}).setOptions({
  title: 'Harmonized Land Surface Temperature Distribution - LST (°C)',
  hAxis: {title: 'Temperature (°C)'},
  vAxis: {title: 'Pixel Count'},
  colors: ['#2b5c8f', '#d95f02']
});
print(chartLST);

// B. NDVI Histogram
var chartNDVI = ui.Chart.image.histogram({
  image: ee.Image.cat([ndvi2014.rename('NDVI 2014'), ndvi2024.rename('NDVI 2024')]),
  region: roi,
  scale: 30,
  maxPixels: 1e9
}).setOptions({
  title: 'Normalized Difference Vegetation Index Distribution - NDVI',
  hAxis: {title: 'NDVI'},
  vAxis: {title: 'Pixel Count'},
  colors: ['#2b5c8f', '#d95f02']
});
print(chartNDVI);

// C. NDBI Histogram
var chartNDBI = ui.Chart.image.histogram({
  image: ee.Image.cat([ndbi2014.rename('NDBI 2014'), ndbi2024.rename('NDBI 2024')]),
  region: roi,
  scale: 30,
  maxPixels: 1e9
}).setOptions({
  title: 'Normalized Difference Built-up Index Distribution - NDBI',
  hAxis: {title: 'NDBI'},
  vAxis: {title: 'Pixel Count'},
  colors: ['#2b5c8f', '#d95f02']
});
print(chartNDBI);

// D. ECI Histogram
var chartECI = ui.Chart.image.histogram({
  image: ee.Image.cat([eci2014.rename('ECI 2014'), eci2024.rename('ECI 2024')]),
  region: roi,
  scale: 30,
  maxPixels: 1e9
}).setOptions({
  title: 'Environmental Criticality Index Distribution - ECI',
  hAxis: {title: 'ECI (0 to 1)'},
  vAxis: {title: 'Pixel Count'},
  colors: ['#2b5c8f', '#d95f02']
});
print(chartECI);

// ==============================================================================
// 7. EXPORTING 8 FILES TO GOOGLE DRIVE (UTM Zone 31N - EPSG:32631)
// ==============================================================================

// 1 & 2. Land Surface Temperature (LST)
Export.image.toDrive({
  image: lst2014, description: 'Parakou_LST_2014', folder: 'GEE_Parakou',
  region: roi, scale: 30, crs: 'EPSG:32631', maxPixels: 1e9
});
Export.image.toDrive({
  image: lst2024, description: 'Parakou_LST_2024', folder: 'GEE_Parakou',
  region: roi, scale: 30, crs: 'EPSG:32631', maxPixels: 1e9
});

// 3 & 4. Vegetation Index (NDVI)
Export.image.toDrive({
  image: ndvi2014, description: 'Parakou_NDVI_2014', folder: 'GEE_Parakou',
  region: roi, scale: 30, crs: 'EPSG:32631', maxPixels: 1e9
});
Export.image.toDrive({
  image: ndvi2024, description: 'Parakou_NDVI_2024', folder: 'GEE_Parakou',
  region: roi, scale: 30, crs: 'EPSG:32631', maxPixels: 1e9
});

// 5 & 6. Built-up Index (NDBI)
Export.image.toDrive({
  image: ndbi2014, description: 'Parakou_NDBI_2014', folder: 'GEE_Parakou',
  region: roi, scale: 30, crs: 'EPSG:32631', maxPixels: 1e9
});
Export.image.toDrive({
  image: ndbi2024, description: 'Parakou_NDBI_2024', folder: 'GEE_Parakou',
  region: roi, scale: 30, crs: 'EPSG:32631', maxPixels: 1e9
});

// 7 & 8. Environmental Criticality Index (ECI)
Export.image.toDrive({
  image: eci2014, description: 'Parakou_ECI_2014', folder: 'GEE_Parakou',
  region: roi, scale: 30, crs: 'EPSG:32631', maxPixels: 1e9
});
Export.image.toDrive({
  image: eci2024, description: 'Parakou_ECI_2024', folder: 'GEE_Parakou',
  region: roi, scale: 30, crs: 'EPSG:32631', maxPixels: 1e9
});