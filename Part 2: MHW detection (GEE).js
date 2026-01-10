/*
EnviroChange: Marine Heatwaves: An overlooked Driver of Baltic Cyanobacteria Blooms?
Part 2: Detecting Marine Heatwaves (MHW) in the Balitc Sea  

This script detects  Marine Heatwave (MHW)conditions based on daily sea surface temperature (SST) data 
from the NOAA OISST v2.1 dataset for summer months (June, July, August) across the Baltic Sea ROI. 
The period covered is 2015–2025.

Worksteps
- Load daily SST data (NOAA OISST v2.1) and convert values to °C.
- Restrict the dataset spatially to the region of interest (ROI) and temporally to summer months.
- Construct a daily SST climatology for the baseline period 1991–2020 using a
  day-of-year (DOY) approach (Smith et al. 2025).
- Compute a 90th percentile SST threshold for each summer DOY (robust to leap years).
- Match each daily SST observation (2015–2025) with its corresponding climatological
  threshold.
- Aggregate spatial mean SST and threshold values over the ROI for each day.
- Export a CSV for further analysis in R (See GitHub-Repository).
  
Key assumptions and notes
- ROI must be provided as a geometry or FeatureCollection (variable: roi).
- Marine heatwaves are identified from surface SST only
- The use of a fixed 1991–2020 baseline follows standard MHW detection (Smith et al. 2025).
- The script does not directly detect events: event duration and intensity metrics
  are calculated in post-processing (R, See Github Repository).
- The exported time series is used for marine heatwave event detection and intensity
  analysis in R following Hobday et al. 2016.

Author: Justin Lingg-Laham and Jonah van den Bos
*/

// Region of Interest (ROI)
var roiFC   = ee.FeatureCollection(roi);
var roiGeom = roiFC.geometry().buffer(25000); // Buffer (25km) to compensate for the low OISST resolution (ca.0.25°). Ensures that pixel centers lie within the geometry (data extraction).

// Configuration (Timescle, Resolution, Export)
var startYear = 2015;                // analysis start (Satellite-Images)
var endYear   = 2025;                // analysis end   (Satellite-Images)
var months    = [6, 7, 8];           // Summerperiod (June, July, August)
var baseStart = 1991;                // start (climatology Baseline like Smith et al. 2025)
var baseEnd   = 2020;                // end (climatology Baseline like Smith et al. 2025)
var scale     = 20000;               // Spatial Resolution (m)
var outFolder = "GEE_Thesis_Data";   // export result as CSV to Google Drive 

// NOAA-OISST ImageCollection
var oisst = ee.ImageCollection("NOAA/CDR/OISST/V2_1")
  .filterBounds(roiGeom)
  .select("sst")
  .map(function(img){
    return img.multiply(0.01) // scaled sst-values by 0.01 to convert to °C
              .copyProperties(img, ["system:time_start"]);
  });

// construct Baseline Period: 1991-2020(Smith et al. 2025).
var baseline = oisst
  .filter(ee.Filter.calendarRange(baseStart, baseEnd, "year"))
  .filter(ee.Filter.calendarRange(6, 8, "month")) // only summer months 
  .map(function(img){
    var doy = img.date().getRelative("day", "year").add(1);
    return img.set("doy", doy);  // using day-of-year approach
  });

// Robust Day-Of-Year range for summer months (accounts for leap years)
var doyList = ee.List.sequence(140, 250);

// Daily 90th percentile SST Threshold for each DOY
var climImages = doyList.map(function(d){
  d = ee.Number(d);
  var col = baseline.filter(ee.Filter.eq("doy", d));
  var p90 = ee.Image(ee.Algorithms.If(
    col.size().gt(0),
    col.reduce(ee.Reducer.percentile([90])).rename("threshold"),
    ee.Image.constant(-9999).rename("threshold")
  ));
  return p90.updateMask(p90.neq(-9999)).set("doy", d); 
});

// Output climImages as ImageCollection (robust)
var climCollection = ee.ImageCollection(ee.List(climImages)).map(function(img){
  return ee.Image(img).set("doy", ee.Image(img).get("doy")); 
});

// Daily SST time series (2015-2025) for summer
var dailyImgs = oisst
  .filter(ee.Filter.calendarRange(startYear, endYear, "year"))
  .filter(ee.Filter.calendarRange(6, 8, "month"));

var dailyFeatures = dailyImgs.map(function(img){ 
  var doy = ee.Number(img.date().getRelative("day", "year")).add(1);
  // each daily SST matches its climatological Treshold
  var matched = climCollection.filter(ee.Filter.eq("doy", doy));
  var thrImg = ee.Image(matched.first());

  // calculate Metrics over ROI
  var stats = ee.Dictionary(ee.Algorithms.If( // only if data exists
    matched.size().gt(0), 
    img.addBands(thrImg).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roiGeom,
      scale: scale,
      bestEffort: true,
      tileScale: 4
    }),
    // if no Data: just SST
    img.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: roiGeom,
      scale: scale,
      bestEffort: true,
      tileScale: 4
    })
  ));

  var sstVal = stats.get("sst");
  var thrVal = stats.get("threshold");

  return ee.Feature(null, {
    date: img.date().format("yyyy-MM-dd"),
    sst: sstVal,
    threshold: thrVal,
    region: "roi"
  });
});

// Export CSV
Export.table.toDrive({
  collection: ee.FeatureCollection(dailyFeatures),
  description: "Daily_SST_p90_JJA_2015_2025_roi",
  folder: outFolder,
  fileFormat: "CSV",
  selectors: ["date", "sst", "threshold", "region"]
});

print("preview features", ee.FeatureCollection(dailyFeatures).limit(10));

// Debug: random sample Date
var sampleDate = ee.Date("2015-06-01");
var sampleImg = oisst.filterDate(sampleDate, sampleDate.advance(1, "day")).first();
var doy = ee.Number(sampleDate.getRelative("day", "year")).add(1);
var matched = climCollection.filter(ee.Filter.eq("doy", doy));
var thrImgRaw = ee.Image(ee.Algorithms.If(matched.size().gt(0), matched.first(), ee.Image.constant(-9999).rename("threshold")));

var thrDirect = thrImgRaw.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: roiGeom,
  scale: scale,
  bestEffort: true,
  tileScale: 1,
  maxPixels: 1e9
});

// reduceRegion with tileScale 4 
var thrTile4 = thrImgRaw.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: roiGeom,
  scale: scale,
  bestEffort: true,
  tileScale: 4,
  maxPixels: 1e9
});

// thrImg.clipToCollection(roiFC) + reduceRegion
var thrClipped = thrImgRaw.clipToCollection(roiFC).reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: roiGeom, 
  scale: scale,
  bestEffort: true,
  tileScale: 1,
  maxPixels: 1e9
});

// thrImg reduced on roiFC geometry (use roiFC.geometry())
var thrGeom = thrImgRaw.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: roiFC.geometry(),
  scale: scale,
  bestEffort: true,
  tileScale: 1,
  maxPixels: 1e9
});

// sst value for comparison
var sstDirect = sampleImg.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: roiGeom,
  scale: scale,
  bestEffort: true,
  tileScale: 1,
  maxPixels: 1e9
});

print("matched.size", matched.size());
print("matched.first (info)", matched.first());
print("thrImgRaw info", thrImgRaw);
print("thrDirect (tileScale=1)", thrDirect);
print("thrTile4 (tileScale=4)", thrTile4);
print("thrClipped (clipToCollection)", thrClipped);
print("thrGeom (reduceRegion with roiFC.geometry())", thrGeom);
print("sstDirect", sstDirect);
