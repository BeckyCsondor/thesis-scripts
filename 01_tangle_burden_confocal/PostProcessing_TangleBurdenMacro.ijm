// ===== CREATE LABELS & ANALYZE TANGLES =====
// Creates labeled map from binary, then analyzes

// --- USER SETTINGS ---
px = 0.577;      // voxel size X µm
py = 0.577;      // voxel size Y µm
pz = 0.4985;     // voxel size Z µm
minVol_um3 = 10; // minimum tangle volume to include
maxVol_um3 = 200; // maximum tangle volume (exclude large outliers)
hippoVol = 5229750; // hippocampus volume from main analysis (µm³)
// -------------------------------------------

print("\\Clear");
print("=================================");
print("LABEL & ANALYZE TANGLES");
print("=================================");

// ====================================================================
// STEP 1: FIND OR CREATE LABELED MAP
// ====================================================================

// Look for existing labeled map
labelMapTitle = "";
binaryMapTitle = "";
imgList = getList("image.titles");

print("Available images:");
for (i = 0; i < imgList.length; i++) {
    print("  " + (i+1) + ". " + imgList[i]);
    
    // Look for labeled map
    if (indexOf(imgList[i], "-lbl") >= 0) {
        labelMapTitle = imgList[i];
    }
    
    // Look for binary tangle map
    if (indexOf(imgList[i], "tauBinary") >= 0 || 
        indexOf(imgList[i], "tangleMap") >= 0) {
        binaryMapTitle = imgList[i];
    }
}

// If we have a labeled map already, use it
if (labelMapTitle != "") {
    print("\n✓ Found labeled map: " + labelMapTitle);
    selectWindow(labelMapTitle);
    
} else if (binaryMapTitle != "") {
    // Create labeled map from binary
    print("\n→ Creating labeled map from: " + binaryMapTitle);
    selectWindow(binaryMapTitle);
    
    // Make sure it's binary
    getStatistics(area, mean, min, max);
    if (max > 255) {
        print("  Converting to 8-bit...");
        run("8-bit");
    }
    
    print("  Running Connected Components Labeling...");
    print("  (This labels each separate object with a unique number)");
    
    run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
    wait(2000);
    
    // The labeled image should now be active
    labelMapTitle = getTitle();
    print("  ✓ Created: " + labelMapTitle);
    
    // Apply color LUT for visualization
    run("Duplicate...", "title=TangleMap_Color duplicate");
    run("glasbey_on_dark");
    run("Enhance Contrast", "saturated=0.35");
    print("  ✓ Color version: TangleMap_Color");
    
    // Go back to the labeled map for analysis
    selectWindow(labelMapTitle);
    
} else {
    print("\nERROR: No tangle binary map found!");
    print("Please run the main analysis first to create tauBinary or tangleMap");
    exit();
}

// ====================================================================
// STEP 2: ANALYZE THE LABELED MAP
// ====================================================================
print("\n=== ANALYZING LABELED REGIONS ===");
selectWindow(labelMapTitle);

print("Running 3D region analysis...");
print("(Measuring volume, shape, etc. for each tangle)");

run("Analyze Regions 3D", "volume surface_area sphericity euler_number centroid");

wait(2000);

// The results are in a table named after the label map
resultsTableName = labelMapTitle + "-morpho";

// Check if it exists
allWindows = getList("window.titles");
tableExists = false;
for (i = 0; i < allWindows.length; i++) {
    if (allWindows[i] == resultsTableName) {
        tableExists = true;
        break;
    }
}

if (!tableExists) {
    print("\nERROR: Analysis results table not found!");
    print("Expected: " + resultsTableName);
    print("\nAvailable windows:");
    for (i = 0; i < allWindows.length; i++) {
        print("  " + allWindows[i]);
    }
    exit();
}

// Select the morpho results table
print("\nReading from: " + resultsTableName);
selectWindow(resultsTableName);

// Use Table API to get row count
n = Table.size(resultsTableName);

print("Total labeled objects: " + n);

if (n == 0) {
    print("\nERROR: No objects detected!");
    print("The labeled map appears to be empty or all background");
    exit();
}

// Arrays to store data
volumes = newArray(n);
sphericities = newArray(n);
surfaces = newArray(n);
validCount = 0;

totalVolume = 0;

for (i = 0; i < n; i++) {
    vol = getResult("Volume", i);
    
    // Check if volume is in voxels (values < 1000 suggest voxel units)
    // Convert to µm³ if needed
    if (!isNaN(vol)) {
        // If volumes are very small, likely in voxels - convert to µm³
        if (vol < 1000) {
            vol = vol * px * py * pz;  // Convert voxels to µm³
        }
        
        if (vol >= minVol_um3 && vol <= maxVol_um3) {
            volumes[validCount] = vol;
            sphericities[validCount] = getResult("Sphericity", i);
            surfaces[validCount] = getResult("SurfaceArea", i);
            
            totalVolume += vol;
            validCount++;
        }
    }
}

if (validCount == 0) {
    print("\nERROR: No tangles passed the size filter!");
    print("Size filter: " + minVol_um3 + " - " + maxVol_um3 + " µm³");
    print("All " + n + " detected objects were outside this range");
    print("Try adjusting minVol_um3 or maxVol_um3 at the top of the script");
    exit();
}

// Trim arrays
volumes = Array.trim(volumes, validCount);
sphericities = Array.trim(sphericities, validCount);
surfaces = Array.trim(surfaces, validCount);

filtered = n - validCount;
if (filtered > 0) {
    print("Filtered out: " + filtered + " objects (outside " + minVol_um3 + "-" + maxVol_um3 + " µm³ range)");
    
    // Count how many were too small vs too large
    tooSmall = 0;
    tooLarge = 0;
    for (i = 0; i < n; i++) {
        vol = getResult("Volume", i);
        if (!isNaN(vol)) {
            if (vol < 1000) vol = vol * px * py * pz;
            if (vol < minVol_um3) tooSmall++;
            else if (vol > maxVol_um3) tooLarge++;
        }
    }
    print("  Too small (<" + minVol_um3 + "): " + tooSmall);
    print("  Too large (>" + maxVol_um3 + "): " + tooLarge);
}
print("Valid tangles for analysis: " + validCount);

// ====================================================================
// STEP 4: CALCULATE STATISTICS
// ====================================================================
print("\n=== VOLUME STATISTICS ===");

Array.sort(volumes);

minVol = volumes[0];
maxVol = volumes[validCount - 1];
medianVol = volumes[floor(validCount / 2)];

p10 = volumes[floor(validCount * 0.10)];
p25 = volumes[floor(validCount * 0.25)];
p75 = volumes[floor(validCount * 0.75)];
p90 = volumes[floor(validCount * 0.90)];
p95 = volumes[floor(validCount * 0.95)];
p99 = volumes[floor(validCount * 0.99)];

avgVol = totalVolume / validCount;

// Standard deviation
sumSqDiff = 0;
for (i = 0; i < validCount; i++) {
    diff = volumes[i] - avgVol;
    sumSqDiff += diff * diff;
}
stdDev = sqrt(sumSqDiff / validCount);

print("Count: " + validCount + " tangles");
print("");
print("Min:        " + d2s(minVol, 2) + " µm³");
print("10th %ile:  " + d2s(p10, 2) + " µm³");
print("25th %ile:  " + d2s(p25, 2) + " µm³");
print("Median:     " + d2s(medianVol, 2) + " µm³");
print("Mean:       " + d2s(avgVol, 2) + " µm³");
print("Std Dev:    " + d2s(stdDev, 2) + " µm³");
print("75th %ile:  " + d2s(p75, 2) + " µm³");
print("90th %ile:  " + d2s(p90, 2) + " µm³");
print("95th %ile:  " + d2s(p95, 2) + " µm³");
print("99th %ile:  " + d2s(p99, 2) + " µm³");
print("Max:        " + d2s(maxVol, 2) + " µm³");

// ====================================================================
// SIZE CATEGORIES
// ====================================================================
print("\n=== SIZE DISTRIBUTION ===");

debris = 0;      // < 15 µm³
small = 0;       // 15-30 µm³
medium = 0;      // 30-60 µm³
large = 0;       // 60-120 µm³
veryLarge = 0;   // > 120 µm³

for (i = 0; i < validCount; i++) {
    v = volumes[i];
    if (v < 15) debris++;
    else if (v < 30) small++;
    else if (v < 60) medium++;
    else if (v < 120) large++;
    else veryLarge++;
}

print("< 15 µm³:    " + debris + " (" + d2s(100*debris/validCount, 1) + "%) - debris?");
print("15-30 µm³:   " + small + " (" + d2s(100*small/validCount, 1) + "%)");
print("30-60 µm³:   " + medium + " (" + d2s(100*medium/validCount, 1) + "%)");
print("60-120 µm³:  " + large + " (" + d2s(100*large/validCount, 1) + "%)");
print("> 120 µm³:   " + veryLarge + " (" + d2s(100*veryLarge/validCount, 1) + "%) - merged?");

// ====================================================================
// SHAPE ANALYSIS
// ====================================================================
print("\n=== SHAPE ANALYSIS (Sphericity) ===");
print("1.0 = perfect sphere, lower = elongated/irregular");

Array.sort(sphericities);
avgSpher = 0;
for (i = 0; i < validCount; i++) {
    avgSpher += sphericities[i];
}
avgSpher = avgSpher / validCount;
medianSpher = sphericities[floor(validCount / 2)];

print("Median sphericity: " + d2s(medianSpher, 3));
print("Mean sphericity:   " + d2s(avgSpher, 3));

spherical = 0;
elongated = 0;
irregular = 0;

for (i = 0; i < validCount; i++) {
    s = sphericities[i];
    if (s > 0.8) spherical++;
    else if (s > 0.5) elongated++;
    else irregular++;
}

print("  Spherical (>0.8):   " + spherical + " (" + d2s(100*spherical/validCount, 1) + "%)");
print("  Elongated (0.5-0.8): " + elongated + " (" + d2s(100*elongated/validCount, 1) + "%)");
print("  Irregular (<0.5):   " + irregular + " (" + d2s(100*irregular/validCount, 1) + "%)");

// ====================================================================
// BURDEN CALCULATIONS
// ====================================================================
print("\n=== BURDEN & DENSITY ===");
print("Hippocampus volume: " + d2s(hippoVol, 0) + " µm³");
print("Total tangle volume: " + d2s(totalVolume, 2) + " µm³");
print("");

burden = (totalVolume / hippoVol) * 100;
density = (validCount / hippoVol) * 1e9;

print("╔═══════════════════════════════════╗");
print("║  BURDEN:  " + d2s(burden, 4) + " %             ║");
print("║  DENSITY: " + d2s(density, 1) + " /mm³      ║");
print("╚═══════════════════════════════════╝");

// ====================================================================
// CREATE VISUALIZATIONS
// ====================================================================
print("\n=== CREATING PLOTS ===");

// HISTOGRAM
histMax = p99 * 1.2;
binWidth = histMax / 50;

bins = newArray(50);
for (i = 0; i < 50; i++) bins[i] = 0;

for (i = 0; i < validCount; i++) {
    v = volumes[i];
    if (v < histMax) {
        binIndex = floor(v / binWidth);
        if (binIndex >= 50) binIndex = 49;
        bins[binIndex]++;
    }
}

Array.getStatistics(bins, binMin, binMax);

Plot.create("Tangle Volume Histogram", "Volume (µm³)", "Count");
Plot.setFontSize(14);

xValues = newArray(50);
for (i = 0; i < 50; i++) {
    xValues[i] = i * binWidth + binWidth/2;
}

Plot.setColor("blue");
Plot.add("bar", xValues, bins);

Plot.setColor("red");
Plot.setLineWidth(3);
Plot.drawLine(medianVol, 0, medianVol, binMax);
Plot.addText("Median: " + d2s(medianVol, 1), medianVol * 1.05, binMax * 0.92);

Plot.setColor("orange");
Plot.setLineWidth(3);
Plot.drawLine(avgVol, 0, avgVol, binMax * 0.85);
Plot.addText("Mean: " + d2s(avgVol, 1), avgVol * 1.05, binMax * 0.77);

Plot.setLimits(0, histMax, 0, binMax * 1.15);
Plot.show();

print("✓ Histogram created");

// BOX PLOT - Fixed for single sample display
Plot.create("Tangle Volume Distribution", "Tangle Size Categories", "Volume (µm³)");
Plot.setFontSize(14);

// Create a proper box plot showing the distribution
boxLeft = 0.8;
boxRight = 1.2;
boxCenter = 1.0;

// Whiskers (10th to 90th percentile)
Plot.setColor("black");
Plot.setLineWidth(2);
Plot.drawLine(boxCenter, p10, boxCenter, p25);
Plot.drawLine(boxCenter, p75, boxCenter, p90);
Plot.drawLine(boxLeft - 0.1, p10, boxRight + 0.1, p10);  // whisker caps
Plot.drawLine(boxLeft - 0.1, p90, boxRight + 0.1, p90);

// Box (25th to 75th)
Plot.setColor("cyan");
xBox = newArray(boxLeft, boxRight, boxRight, boxLeft, boxLeft);
yBox = newArray(p25, p25, p75, p75, p25);
Plot.add("filled", xBox, yBox);

Plot.setColor("black");
Plot.setLineWidth(2);
Plot.add("line", xBox, yBox);

// Median line (red)
Plot.setColor("red");
Plot.setLineWidth(4);
Plot.drawLine(boxLeft, medianVol, boxRight, medianVol);

// Mean dot (orange)
Plot.setColor("orange");
Plot.setLineWidth(3);
xMean = newArray(1);
yMean = newArray(1);
xMean[0] = boxCenter;
yMean[0] = avgVol;
Plot.add("circle", xMean, yMean);

// Add text labels
Plot.setColor("black");
Plot.addText("10%: " + d2s(p10, 1), 1.3, p10);
Plot.addText("25%: " + d2s(p25, 1), 1.3, p25);
Plot.addText("Median: " + d2s(medianVol, 1), 1.3, medianVol);
Plot.addText("Mean: " + d2s(avgVol, 1), 1.3, avgVol);
Plot.addText("75%: " + d2s(p75, 1), 1.3, p75);
Plot.addText("90%: " + d2s(p90, 1), 1.3, p90);

Plot.setLimits(0.5, 2, 0, maxVol * 1.1);
Plot.show();

print("✓ Box plot created");

// SCATTER
Plot.create("Volume vs Sphericity", "Volume (µm³)", "Sphericity");
Plot.setFontSize(14);
Plot.setColor("blue");
Plot.add("circle", volumes, sphericities);
Plot.setLimits(0, p99 * 1.1, 0, 1.1);
Plot.show();

print("✓ Scatter plot created");

// ====================================================================
// SUMMARY
// ====================================================================
print("\n=================================");
print("✓ ANALYSIS COMPLETE!");
print("=================================");
print("\nMedian tangle: " + d2s(medianVol, 1) + " µm³");
print("Mean tangle: " + d2s(avgVol, 1) + " µm³");
print("Burden: " + d2s(burden, 3) + " %");
print("Density: " + d2s(density, 0) + " tangles/mm³");
print("=================================");

// ====================================================================
// SAVE TO SUMMARY TABLE (for multiple samples)
// ====================================================================
print("\n=== SAVING TO SUMMARY TABLE ===");

// Get sample name from user or use image title
Dialog.create("Sample Name");
Dialog.addString("Sample/Image name:", labelMapTitle);
Dialog.addMessage("This will be saved to a summary table");
Dialog.show();
sampleName = Dialog.getString();

// Create or update summary table
if (!isOpen("Tangle_Analysis_Summary")) {
    // Create new table
    Table.create("Tangle_Analysis_Summary");
    Table.setColumn("Sample", newArray(0));
}

// Add this sample's data
selectWindow("Tangle_Analysis_Summary");
row = Table.size;

Table.set("Sample", row, sampleName);
Table.set("Hippo_Vol_um3", row, hippoVol);
Table.set("Tangle_Count", row, validCount);
Table.set("Total_Tangle_Vol_um3", row, totalVolume);
Table.set("Burden_%", row, burden);
Table.set("Density_per_mm3", row, density);
Table.set("Median_Vol_um3", row, medianVol);
Table.set("Mean_Vol_um3", row, avgVol);
Table.set("StdDev_um3", row, stdDev);
Table.set("Mean_Sphericity", row, avgSpher);
Table.set("Debris_<15um3", row, debris);
Table.set("Small_15-30", row, small);
Table.set("Medium_30-60", row, medium);
Table.set("Large_60-120", row, large);
Table.set("VeryLarge_>120", row, veryLarge);

Table.update;

print("✓ Added to summary table: " + sampleName);
print("\nYou can:");
print("  • Run this script again for more samples");
print("  • Summary table will accumulate all results");
print("  • File → Save As... to export the summary table");
print("=================================");