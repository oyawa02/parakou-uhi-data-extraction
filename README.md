# Urban Heat Islands Mapping in Parakou (Benin) using Satellite Remote Sensing

Google Earth Engine (GEE) JavaScript scripts developed to evaluate urban microclimate changes, urban heat islands (UHI), and environmental criticality between 2014 and 2024 in Parakou, Benin.

## Overview
This repository contains the complete workflow for:
1. **Preprocessing Landsat 8 Collection 2 Level-2 imagery** (Cloud masking, scaling factors).
2. **Calculating key spectral indices**: 
   - LST (Land Surface Temperature in °C) with interannual thermal bias correction.
   - NDVI (Normalized Difference Vegetation Index).
   - NDBI (Normalized Difference Built-up Index).
   - ECI (Environmental Criticality Index).
3. **Generating comparative charts** and exporting 8 raster datasets to Google Drive (UTM Zone 31N - EPSG:32631).

## Prerequisites
- A free **Google Earth Engine (GEE)** account ([Sign up here](https://earthengine.google.com/)).
- Access to an administrative boundary asset for Parakou (update the `roi` asset path in the script if you use your own vector asset).

## Usage
1. Open the [Google Earth Engine Code Editor](https://code.earthengine.google.com/).
2. Create a new script and paste the contents of `parakou_uhi_analysis.js`.
3. Run the script to visualize the layers and charts.
4. Go to the **Tasks** tab in the GEE interface to run the 8 exports to your Google Drive.

## Authors
- **Ozias OROU YAWA** (Population Dynamics Lab POP Dyn, University of Parakou)
- **Tanguy Mesreil HOUEGBE** (BiPLAC, LEB, University of Parakou)

## License
This project is licensed under the [MIT License](LICENSE).
