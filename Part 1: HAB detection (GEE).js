/*
EnviroChange: Marine Heatwaves: An overlooked Driver of Baltic Cyanobacteria Blooms?
Part 1: Detecting Phycocyanin in the Baltic Sea 

This script calculates a Phycocyanin (PC) proxy index from Sentinel‑2 Level‑2 (SR_HARMONIZED)
for summer months (June, July, August) across the Baltic Sea ROI. The period covered is 2015–2025.

Worksteps
- Load Satellite-Data (Sentinel‑2 SR_HARMONIZED) and applie a QA + scene classification mask to retain water pixels
  and remove clouds/shadows (Wright et al. 2024).
- Compute a PC proxy per image from spectral bands (B4, B5, B8A) using Methods from Nesheli et al. 2024
- Aggregate results monthly (Summerperiod) for each year in the interval:
  - image_count (number of used S2 images)
  - PC_mean, PC_median, PC_max (spatial average over the ROI of the monthly mean/median/max)
  - PC_p90 (90th percentile of PC)
  - bloom_fraction: fraction of images/pixels in the month with PC > Threshold (Nesheli et al. 2024)
  - bloom_count_approx: bloom_fraction × image_count (approximate bloom event count)
- Export a CSV for further analysis in R (See GitHub-Repository).

Key assumptions and notes
- ROI must be provided as a geometry or FeatureCollection (variable: roi).
- The PC index is a spectral proxy (not a direct in situ measurement); calibration/validation against field data is recommended.
- A fixed bloom threshold (threshold, µg/L) (Nesheli et al. 2024) is used as a specific calibration.

Author: Justin Lingg-Laham and Jonah van den Bos 
*/

// Region of Interest (ROI)
var roiFC = ee.FeatureCollection(roi);
var roiGeom = roiFC.geometry();

// Configuartion (Timescale, Cloud-Treshold, Resolution)
var startYear = 2015;              // Start (satellite-images)
var endYear = 2025;                // End   (satellite-Images)
var months = [6, 7, 8];            // Summerperiod
var s2CloudThreshold = 30;         // Cloud-Treshold (not too low to get enough images)
var reduceScale = 1000;            // Resolution
var tileScale = 4;                 // tileScale = 4 = less storage Space for GEE
var threshold = 10;                // µg/L bloom threshold (moderat) (Nesheli et al. 2024, table 4)
var testMode = false;              // Debug

var exportCSV = true;              // export result as CSV to Google Drive 
var exportFolder = "GEE_Thesis_Data";

// Function: CloudMask (Wright et al. 2024)
function maskS2clouds(image) {
  var qa = image.select("QA60");
  var scl = image.select("SCL");
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var qaMask = qa.bitwiseAnd(cloudBitMask).eq(0)
               .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  var waterMask = scl.eq(6);
  var cloudMask = scl.eq(8).or(scl.eq(9)).or(scl.eq(10));
  var shadowMask = scl.eq(3);
  return image.updateMask(qaMask.and(waterMask).and(cloudMask.not()).and(shadowMask.not()));
}

// Function: Phycocyanin-Detection (Nesheli et al. 2024)
function computePC(image) {
  var scaled = image.select(["B2","B3","B4","B5","B8A"]).multiply(0.0001);
  var R4 = scaled.select("B4"), R5 = scaled.select("B5"), R8A = scaled.select("B8A");
  var baseline705 = R4.add(R8A.subtract(R4).multiply(0.2));
  var MPH = R5.subtract(baseline705).rename("MPH");
  var PC = MPH.multiply(3047.6).add(0.626).max(0).rename("PC_ugL"); //One-Day-Modell for Phycocyanin (Nesheli et al. 2024, table 8)
  return image.addBands([PC]).copyProperties(image, ["system:time_start"]);
}

// Sentinel-2 ImageCollection 
var s2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
  .filterBounds(roiGeom)
  .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", s2CloudThreshold))
  .map(maskS2clouds)
  .map(computePC)
  .select(["PC_ugL"]);

print("S2 base count (server-side):", s2.size());

// Debug: reduce calculation time (with var testMode)
var yearList = ee.List.sequence(startYear, endYear);
if (testMode) { yearList = yearList.slice(0, 2); } 

// nested loop: Data for every Month
var monthPairs = yearList.map(function(y){ // Outer-Loop: go trough every Year
  y = ee.Number(y);
  return ee.List(months).map(function(m){
    m = ee.Number(m);
    var start = ee.Date.fromYMD(y, m, 1);  // Inner-Loop: Go for every of Outer-Loop-Years trough Months
    var end = start.advance(1, "month");

    // Filter: Sentinel-2 Images for that Year and Months + ROI
    var s2month = s2.filterDate(start, end).filterBounds(roiGeom);
    var imageCount = s2month.size();
    var hasImages = imageCount.gt(0); // if hasImages -> compute Metrics (Debug)

    // PC-Mean
    var pc_mean = ee.Algorithms.If(hasImages,
      s2month.mean().select("PC_ugL").reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: roiGeom,
        scale: reduceScale,
        tileScale: tileScale,
        bestEffort: true,
        maxPixels: 1e9
      }).get("PC_ugL"),
      null);
    // PC-Median
    var pc_median = ee.Algorithms.If(hasImages,
      s2month.median().select("PC_ugL").reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: roiGeom,
        scale: reduceScale,
        tileScale: tileScale,
        bestEffort: true,
        maxPixels: 1e9
      }).get("PC_ugL"),
      null);
    //PC-Max
    var pc_max = ee.Algorithms.If(hasImages,
      s2month.max().select("PC_ugL").reduceRegion({
        reducer: ee.Reducer.max(),
        geometry: roiGeom,
        scale: reduceScale,
        tileScale: tileScale,
        bestEffort: true,
        maxPixels: 1e9
      }).get("PC_ugL"),
      null);

    // 90. Percentile
    var pc_p90 = ee.Algorithms.If(hasImages,
      ee.Image(s2month.select("PC_ugL").reduce(ee.Reducer.percentile([90]))).reduceRegion({
        reducer: ee.Reducer.mean(),
        geometry: roiGeom,
        scale: reduceScale,
        tileScale: tileScale,
        bestEffort: true,
        maxPixels: 1e9
      }).get("PC_ugL_p90"),
      null);

    // Bloom fraction: map binary images (PC > Treshold) -> mean over time -> reduceRegion(mean)
    var bloomFraction = ee.Algorithms.If(hasImages,
      s2month.map(function(img){ return img.select("PC_ugL").gt(threshold).rename("b"); })
             .select("b").mean()
             .reduceRegion({
               reducer: ee.Reducer.mean(),
               geometry: roiGeom,
               scale: reduceScale,
               tileScale: tileScale,
               bestEffort: true,
               maxPixels: 1e9
             }).get("b"),
      null);

    // Bloom count approx (estimated value) (fraction * imageCount) 
    var bloomCountApprox = ee.Algorithms.If(hasImages,
      ee.Number(bloomFraction).multiply(imageCount),
      0);

    // Feature: All For One
    return ee.Feature(null, {
      "year": y,
      "month": m,
      "start": start.format("YYYY-MM-dd"),
      "end": end.advance(-1, "day").format("YYYY-MM-dd"),
      "image_count": imageCount,
      "PC_mean": pc_mean,
      "PC_median": pc_median,
      "PC_max": pc_max,
      "PC_p90": pc_p90,
      "threshold_used": threshold,
      "bloom_fraction": bloomFraction,
      "bloom_count_approx": bloomCountApprox
    });
  });
}).flatten(); // One Flat List (better for Export as CSV)

// FeatureCollection: MonthPairs
var monthlyFC = ee.FeatureCollection(monthPairs);

// Export CSV
if (exportCSV) {
  Export.table.toDrive({
    collection: monthlyFC,
    description: "S2_PC_monthly_2015_2025_roi",
    folder: exportFolder,
    fileFormat: "CSV"
  });
}
// Check before Export
print("starting check...");

// visual Check
var debugYear = 2020;
var debugMonth = 8; 
var debugStart = ee.Date.fromYMD(debugYear, debugMonth, 1);
var debugEnd = ee.Date.fromYMD(debugYear, 8, 30);
var s2Debug = s2.filterDate(debugStart, debugEnd);

// Number of Images 
s2Debug.size().evaluate(function(count) {
  print("Number S2-Images Summer" + debugYear + "(after Cloudfilter):", count);
});

var visualization = {
  min: 0, 
  max: 30, 
  palette: ["blue", "cyan", "green", "yellow", "red"]
};

Map.centerObject(roiGeom, 10);
Map.addLayer(roiGeom, {color: "red"}, 'ROI', false);
Map.addLayer(s2Debug.mean().clip(roiGeom), visualization, "PC Mean Summer" + debugYear);


// Isolated Datacheck
var start = ee.Date.fromYMD(debugYear, debugMonth, 1);
var end = start.advance(1, "month");
var s2month = s2.filterDate(start, end).filterBounds(roiGeom);
var imageCount = s2month.size();
var hasImages = imageCount.gt(0);

var combinedReducer = ee.Reducer.mean()
  .combine({reducer2: ee.Reducer.median(), sharedInputs: true})
  .combine({reducer2: ee.Reducer.max(), sharedInputs: true});

var timeReducer = ee.Reducer.mean()
  .combine({reducer2: ee.Reducer.median(), sharedInputs: true})
  .combine({reducer2: ee.Reducer.max(), sharedInputs: true});

var reducedImage = ee.Algorithms.If(hasImages,
    s2month.select("PC_ugL").reduce(timeReducer),
    ee.Image(0).rename(["PC_ugL_mean", "PC_ugL_median", "PC_ugL_max"]) 
);

var stats = ee.Dictionary(ee.Algorithms.If(hasImages,
  ee.Image(reducedImage).reduceRegion({ 
    reducer: ee.Reducer.first(), 
    geometry: roiGeom,
    scale: reduceScale,
    tileScale: tileScale,
    bestEffort: true,
    maxPixels: 1e9
  }),
  {"PC_ugL_mean": null, "PC_ugL_median": null, "PC_ugL_max": null}
));

var bloomFraction = ee.Algorithms.If(hasImages,
  s2month.map(function(img){ return img.gt(threshold); })
         .mean() 
         .reduceRegion({
           reducer: ee.Reducer.mean(), 
           geometry: roiGeom,
           scale: reduceScale,
           tileScale: tileScale,
           bestEffort: true,
           maxPixels: 1e9
         }).get("PC_ugL"), 
  null);

// Feature for Test-Month
var testFeature = ee.Feature(null, {
  "year": debugYear,
  "month": debugMonth,
  "PC_mean": stats.get("PC_ugL_mean"),
  "bloom_fraction": bloomFraction
});

print("calculate test-Datapoints for " + debugYear + "/" + debugMonth + ":");
testFeature.evaluate(function(feature) {
  if (feature && feature.properties.PC_mean !== null) {
    print("Sucsess! Test-Datapoint:", feature.properties);
    print("Result is correct! start export");
  } else {
    print("Attention!: The Value is null.");
  }
});
