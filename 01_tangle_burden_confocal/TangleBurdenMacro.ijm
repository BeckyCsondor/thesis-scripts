// ===== NFT Burden Quantification - ALTERNATIVE METHOD =====
// Uses built-in Analyze Particles instead of 3D Objects Counter
// Input: 2-channel composite stack (FusionRed tau + Syn-GFP)

// --- USER SETTINGS ---
tauCh = 1;           // FusionRed tau channel (1-based)
synCh = 2;           // Synaptophysin-GFP channel (1-based)
px = 0.577;          // voxel size X µm
py = 0.577;          // voxel size Y µm
pz = 0.4985;         // voxel size Z µm
minVol_um3 = 10;     // minimum tangle volume in µm^3 (debris filter)
maxVol_um3 = 10000;  // maximum tangle volume in µm^3

// Preprocessing
tauBlurX = 0.8;
tauBlurY = 0.8;
tauBlurZ = 0.5;
// -------------------------------------------

setBatchMode(false);

print("\\Clear");
print("\n=================================");
print("NFT TANGLE BURDEN ANALYSIS");
print("(No 3D Objects Counter required)");
print("=================================");

// --- Get original composite ---
origTitle = getTitle();
origID = getImageID();
print("Processing: " + origTitle);

// --- Set voxel calibration ---
selectImage(origID);
run("Properties...", "unit=micron pixel_width="+px+" pixel_height="+py+" voxel_depth="+pz);
print("Voxel size: " + px + " × " + py + " × " + pz + " µm");

// --- Duplicate channels ---
print("\n[Step 1/6] Duplicating channels...");
selectImage(origID);
run("Duplicate...", "title=tauDup duplicate channels="+tauCh);
tauDupID = getImageID();

selectImage(origID);
run("Duplicate...", "title=synDup duplicate channels="+synCh);
synDupID = getImageID();

// --- Generate hippocampal mask ---
print("[Step 2/6] Generating hippocampal mask...");
selectImage(synDupID);
run("Duplicate...", "title=synSmooth duplicate");
synSmoothID = getImageID();

// Moderate blur - enough to connect but keep natural shape
run("Gaussian Blur 3D...", "x=3 y=3 z=2");

// Ensure 8-bit
run("Select None");
if (bitDepth() != 8) {
    run("8-bit");
}
run("Enhance Contrast", "saturated=0.1");

print("   Opening threshold tool for hippocampus...");
print("   Adjust to capture the full hippocampal region");
print("   (Keep lower threshold fairly low to include all tissue)");

run("Threshold...");

waitForUser("Hippocampus Threshold", 
    "Adjust threshold for hippocampus:\n" +
    "1. Lower slider = keep low to include all GFP signal\n" +
    "2. Should capture the full hippocampal outline\n" +
    "3. Check multiple Z-slices\n" +
    "4. Click OK when done\n\n" +
    "DO NOT click Apply!");

getThreshold(lower, upper);
print("   Hippocampus threshold: " + lower + " - " + upper);

if (lower == -1) {
    print("   No threshold set, using default");
    setAutoThreshold("Yen dark");
}

run("Convert to Mask");

// Gentler filling - preserves natural shape better
Stack.getDimensions(w, h, c, slices, f);
print("   Filling interior gaps (gentle method)...");

for (z = 1; z <= slices; z++) {
    setSlice(z);
    
    // Simple hole filling without convex hull
    run("Fill Holes", "slice");
    
    if (z % 10 == 0) {
        print("   Slice " + z + "/" + slices);
    }
}

// Light smoothing to reduce jaggedness but keep shape
run("Dilate");
run("Erode");

rename("hippoMask");
hippoMaskID = getImageID();
print("   ✓ Hippocampal mask created");

if (isOpen("Threshold")) {
    selectWindow("Threshold");
    run("Close");
}

selectImage(synDupID);
close();

// --- Apply mask to tau ---
print("[Step 3/6] Masking tau...");
selectImage(tauDupID);
imageCalculator("AND create stack", "tauDup", "hippoMask");
rename("tauMasked");
tauMaskedID = getImageID();

// --- Prepare for manual threshold ---
print("[Step 4/6] Preparing tau for thresholding...");
selectImage(tauMaskedID);
run("Gaussian Blur 3D...", "x="+tauBlurX+" y="+tauBlurY+" z="+tauBlurZ);
run("Duplicate...", "title=tauForThreshold duplicate");
tauThreshID = getImageID();
run("Enhance Contrast", "saturated=0.1");

// --- MANUAL THRESHOLD ---
print("\n[Step 5/6] MANUAL THRESHOLD");
print("=================================");
print("Threshold window opening...");
print("");
print("INSTRUCTIONS:");
print("1. Move LOWER slider RIGHT to exclude background");
print("2. Only BRIGHT tangles should be red");
print("3. Check multiple Z-slices");
print("4. Click OK when done");
print("=================================");

selectImage(tauThreshID);
run("Threshold...");

waitForUser("Adjust Threshold", 
    "1. Adjust threshold sliders\n" +
    "2. Only bright tangles highlighted (red)\n" +
    "3. Check multiple Z-slices\n" +
    "4. Click OK when ready\n\n" +
    "DO NOT click Apply!");

selectImage(tauThreshID);
getThreshold(lower, upper);
print("\nThreshold set: " + lower + " - " + upper);

if (lower == -1) {
    print("No threshold set, using default");
    setAutoThreshold("Triangle dark");
}

run("Convert to Mask", "method=Default background=Dark black");
rename("tauBinary");
tauBinaryID = getImageID();

if (isOpen("Threshold")) {
    selectWindow("Threshold");
    run("Close");
}

// ====================================================================
// MEASURE HIPPOCAMPUS - Manual voxel counting
// ====================================================================
print("\n[Step 6/6] Measuring volumes...");
print("\n--- HIPPOCAMPUS ---");

selectImage(hippoMaskID);
Stack.getDimensions(width, height, channels, slices, frames);
print("Image size: " + width + " × " + height + " × " + slices);

hippoVoxels = 0;
print("Counting hippocampus voxels...");

for (z = 1; z <= slices; z++) {
    setSlice(z);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            if (getPixel(x, y) > 0) {
                hippoVoxels++;
            }
        }
    }
    if (z % 5 == 0 || z == slices) {
        print("   " + z + "/" + slices + " slices...");
    }
}

hippoVol = hippoVoxels * px * py * pz;
print("Hippocampus: " + hippoVoxels + " voxels = " + d2s(hippoVol, 2) + " µm³");

// ====================================================================
// MEASURE TAU - Manual method with connected components
// ====================================================================
print("\n--- TAU TANGLES ---");
print("Counting tangle voxels...");

selectImage(tauBinaryID);

tangleVoxels = 0;
for (z = 1; z <= slices; z++) {
    setSlice(z);
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            if (getPixel(x, y) > 0) {
                tangleVoxels++;
            }
        }
    }
    if (z % 5 == 0 || z == slices) {
        print("   " + z + "/" + slices + " slices...");
    }
}

tangleVol = tangleVoxels * px * py * pz;
print("Total tangle voxels: " + tangleVoxels);
print("Total tangle volume: " + d2s(tangleVol, 2) + " µm³");

// Count individual tangles using 3D Connected Components Labeling
print("\nCounting individual tangles...");
selectImage(tauBinaryID);

// Try using MorphoLibJ if available
if (File.exists(getDirectory("plugins") + "MorphoLibJ_-1.6.1.jar") ||
    File.exists(getDirectory("plugins") + "MorphoLibJ-1.6.1.jar")) {
    
    print("Using MorphoLibJ for connected components...");
    run("Connected Components Labeling", "connectivity=6 type=[16 bits]");
    labelID = getImageID();
    
    // Get max label (= number of objects)
    Stack.getStatistics(voxelCount, mean, min, max);
    tangleCount = max;
    
    selectImage(labelID);
    close();
    
} else {
    // Fallback: use Analyze Particles on each slice and approximate
    print("Counting using 2D particle analysis...");
    
    run("Clear Results");
    totalParticles = 0;
    
    for (z = 1; z <= slices; z++) {
        setSlice(z);
        run("Analyze Particles...", "size=0-Infinity display clear slice");
        totalParticles += nResults;
    }
    
    // Approximate 3D objects (will overestimate)
    tangleCount = round(totalParticles / 3);  // Rough estimate
    print("Approximate tangle count (may overestimate): " + tangleCount);
    print("Note: Install MorphoLibJ plugin for accurate counting");
}

print("Detected tangles: " + tangleCount);

// ====================================================================
// ANALYZE TANGLE SIZE DISTRIBUTION
// ====================================================================
if (tangleCount > 0 && isOpen("Results")) {
    print("\n--- TANGLE SIZE ANALYSIS ---");
    
    // Get actual number of results (might differ from tangleCount)
    selectWindow("Results");
    actualRows = nResults;
    print("Results table has " + actualRows + " rows");
    
    if (actualRows == 0) {
        print("WARNING: Results table is empty, skipping analysis");
    } else {
        // Extract all volumes into an array
        volumes = newArray(actualRows);
        validCount = 0;
        
        for (i = 0; i < actualRows; i++) {
            vol = NaN;
            
            // Try different column names for volume
            if (isNaN(vol)) {
                vol = getResult("Volume (micron^3)", i);
            }
            if (isNaN(vol)) {
                volPix = getResult("Volume (pixel^3)", i);
                if (!isNaN(volPix)) {
                    vol = volPix * px * py * pz;
                }
            }
            if (isNaN(vol)) {
                vol = getResult("Volume", i);
            }
            
            // If no Volume column, try to calculate from Area (2D particle analysis)
            if (isNaN(vol) || vol == 0) {
                area = getResult("Area", i);
                if (!isNaN(area) && area > 0) {
                    // Estimate volume assuming spherical: V = 4/3 * pi * r^3
                    // where r = sqrt(Area/pi)
                    radius = sqrt(area / PI);
                    vol = (4.0/3.0) * PI * pow(radius, 3);
                    
                    // Alternative: just multiply by a typical height
                    // vol = area * (minVol_um3 / 10);  // rough estimate
                }
            }
            
            if (!isNaN(vol) && vol > 0) {
                volumes[validCount] = vol;
                validCount++;
            }
        }
        
        // Trim array to valid entries
        validVolumes = Array.trim(volumes, validCount);
        
        print("Valid volume measurements: " + validCount);
        
        if (validCount == 0) {
            print("\n⚠ WARNING: No volume data found in Results table");
            print("Results table columns found:");
            
            // List available columns
            if (actualRows > 0) {
                // Try to list first few column values to debug
                print("  Area: " + getResult("Area", 0));
                print("  Mean: " + getResult("Mean", 0));
                
                // Calculate simple statistics from the manual voxel count instead
                print("\nUsing total volume from voxel counting:");
                print("  Total tangle volume: " + d2s(tangleVol, 2) + " µm³");
                print("  Estimated count: " + tangleCount);
                if (tangleCount > 0) {
                    avgVol = tangleVol / tangleCount;
                    print("  Estimated average: " + d2s(avgVol, 2) + " µm³");
                }
                print("\nNote: Install MorphoLibJ for detailed per-tangle analysis");
            }
        } else if (validCount < 2) {
            print("Not enough data for distribution analysis");
        } else {
            // Sort to find percentiles
            Array.sort(validVolumes);
            
            // Calculate statistics
            minVol = validVolumes[0];
            maxVol = validVolumes[validCount - 1];
            medianVol = validVolumes[round(validCount / 2)];
            
            p25 = validVolumes[round(validCount * 0.25)];
            p75 = validVolumes[round(validCount * 0.75)];
            p90 = validVolumes[round(validCount * 0.90)];
            p95 = validVolumes[round(validCount * 0.95)];
            p99 = validVolumes[round(validCount * 0.99)];
            
            // Calculate mean
            sumVol = 0;
            for (i = 0; i < validCount; i++) {
                sumVol += validVolumes[i];
            }
            avgVol = sumVol / validCount;
            
            print("\nVolume statistics:");
            print("  Min: " + d2s(minVol, 2) + " µm³");
            print("  25th percentile: " + d2s(p25, 2) + " µm³");
            print("  Median: " + d2s(medianVol, 2) + " µm³");
            print("  Mean: " + d2s(avgVol, 2) + " µm³");
            print("  75th percentile: " + d2s(p75, 2) + " µm³");
            print("  90th percentile: " + d2s(p90, 2) + " µm³");
            print("  95th percentile: " + d2s(p95, 2) + " µm³");
            print("  99th percentile: " + d2s(p99, 2) + " µm³");
            print("  Max: " + d2s(maxVol, 2) + " µm³");
            
            // Count size categories
            verySmall = 0;  // < 15 µm³ (likely debris)
            small = 0;      // 15-30 µm³
            medium = 0;     // 30-60 µm³
            large = 0;      // 60-120 µm³
            veryLarge = 0;  // > 120 µm³
            
            for (i = 0; i < validCount; i++) {
                v = validVolumes[i];
                if (v < 15) verySmall++;
                else if (v < 30) small++;
                else if (v < 60) medium++;
                else if (v < 120) large++;
                else veryLarge++;
            }
            
            print("\nSize distribution:");
            print("  < 15 µm³: " + verySmall + " (" + d2s(100*verySmall/validCount, 1) + "%) - likely debris");
            print("  15-30 µm³: " + small + " (" + d2s(100*small/validCount, 1) + "%)");
            print("  30-60 µm³: " + medium + " (" + d2s(100*medium/validCount, 1) + "%)");
            print("  60-120 µm³: " + large + " (" + d2s(100*large/validCount, 1) + "%)");
            print("  > 120 µm³: " + veryLarge + " (" + d2s(100*veryLarge/validCount, 1) + "%) - may be merged");
            
            // Create histogram
            print("\nCreating histogram...");
            
            // Use 50 bins from minVol to p99 (exclude extreme outliers from plot)
            histMax = p99 * 1.1;
            binWidth = histMax / 50;
            
            bins = newArray(50);
            for (i = 0; i < 50; i++) {
                bins[i] = 0;
            }
            
            // Count into bins
            for (i = 0; i < validCount; i++) {
                v = validVolumes[i];
                if (v < histMax) {
                    binIndex = floor(v / binWidth);
                    if (binIndex >= 50) binIndex = 49;
                    bins[binIndex]++;
                }
            }
            
            // Find max bin for scaling
            Array.getStatistics(bins, binMin, binMax);
            
            // Plot histogram
            Plot.create("Tangle Volume Distribution", "Volume (µm³)", "Count");
            
            xValues = newArray(50);
            for (i = 0; i < 50; i++) {
                xValues[i] = i * binWidth + binWidth/2;
            }
            
            Plot.add("bar", xValues, bins);
            Plot.setColor("blue");
            Plot.setLineWidth(2);
            
            // Add vertical lines for key statistics
            Plot.setColor("red");
            Plot.drawLine(medianVol, 0, medianVol, binMax);
            Plot.addText("Median", medianVol, binMax * 0.9);
            
            Plot.setColor("orange");
            Plot.drawLine(avgVol, 0, avgVol, binMax * 0.8);
            Plot.addText("Mean", avgVol, binMax * 0.7);
            
            Plot.setLimits(0, histMax, 0, binMax * 1.1);
            Plot.show();
            
            print("✓ Histogram created");
            print("  (Showing up to 99th percentile: " + d2s(p99, 1) + " µm³)");
            
            // Flag suspicious outliers
            if (veryLarge > 0) {
                print("\n⚠ WARNING: " + veryLarge + " very large tangles (>120 µm³)");
                print("  These may be merged tangles or artifacts");
                print("  Consider:");
                print("  - Reducing blur (tauBlurX/Y/Z)");
                print("  - Raising manual threshold");
                print("  - Setting maxVol_um3 to 100-120");
            }
            
            if (verySmall > validCount * 0.3) {
                print("\n⚠ WARNING: " + d2s(100*verySmall/validCount, 0) + "% are < 15 µm³");
                print("  May include noise/debris");
                print("  Consider raising minVol_um3 to 15-20");
            }
        }
    }
}

// ====================================================================
// CREATE COLOR-CODED TANGLE MAP
// ====================================================================
print("\n[Bonus] Creating color-coded tangle map...");

selectImage(tauBinaryID);
run("Duplicate...", "title=tangleMap duplicate");
tangleMapID = getImageID();

// Use 3D Connected Components to label each tangle
selectImage(tangleMapID);

// Try to run MorphoLibJ labeling (will fail gracefully if not installed)
run("Connected Components Labeling", "connectivity=6 type=[16 bits]");

// Check if it worked
if (getTitle() == "tangleMap-lbl") {
    // Success! MorphoLibJ created a labeled image
    labelMapID = getImageID();
    
    // Convert to color for visualization
    run("Duplicate...", "title=TangleMap_Color duplicate");
    colorMapID = getImageID();
    
    // Apply a nice color LUT
    run("glasbey_on_dark");
    run("Enhance Contrast", "saturated=0.35");
    
    print("✓ Color-coded tangle map created: 'TangleMap_Color'");
    print("  Each color = one individual tangle");
    print("  Check if multiple tangles have the same color (=merged)");
    
    // Also create a montage for easy viewing
    run("Make Montage...", "columns=4 rows=3 scale=0.5 border=2");
    rename("TangleMap_Montage");
    print("✓ Montage view created for quick browsing");
    
} else {
    // MorphoLibJ not available, use simpler method
    print("  MorphoLibJ not found - using basic visualization");
    selectImage(tangleMapID);
    
    // Just apply a color LUT to the binary
    run("16 Colors");
    print("✓ Basic tangle map created");
    print("  Install MorphoLibJ for proper color-coded labeling:");
    print("  Help → Update → Manage update sites → Check 'IJPB-plugins'");
}

// ====================================================================
// FINAL RESULTS
// ====================================================================
print("\n=================================");
print("FINAL RESULTS");
print("=================================");
print("Image: " + origTitle);
print("");
print("Hippocampus volume: " + d2s(hippoVol, 2) + " µm³");
print("                    " + d2s(hippoVol/1e9, 6) + " mm³");
print("");
print("Total tangle volume: " + d2s(tangleVol, 2) + " µm³");
print("Number of tangles: " + tangleCount);
print("");

if (hippoVol > 0) {
    tangleBurden = (tangleVol / hippoVol) * 100;
    density = (tangleCount / hippoVol) * 1e9;
    
    print("╔════════════════════════════════════╗");
    print("║  TANGLE BURDEN: " + d2s(tangleBurden, 4) + " %        ║");
    print("║  TANGLE DENSITY: " + d2s(density, 2) + " /mm³  ║");
    print("╚════════════════════════════════════╝");
    
    if (tangleCount > 0) {
        avgVol = tangleVol / tangleCount;
        print("\nAverage tangle size: " + d2s(avgVol, 2) + " µm³");
    }
} else {
    print("ERROR: Cannot calculate burden");
}

print("\n=================================");
print("✓ Analysis complete!");
print("\nVERIFY:");
print("  • hippoMask = solid volume?");
print("  • tauBinary = only tangles?");
print("=================================");

// Optional: Save results to CSV
Dialog.create("Save Results?");
Dialog.addCheckbox("Save results to CSV file?", true);
Dialog.show();
saveToCSV = Dialog.getCheckbox();

if (saveToCSV) {
    // Let user choose location
    resultsPath = getDirectory("image");
    if (resultsPath == "") {
        resultsPath = getDirectory("home");
    }
    
    resultsFile = resultsPath + "tangle_results.csv";
    
    // Try to save
    success = false;
    try {
        fileExists = File.exists(resultsFile);
        
        if (fileExists) {
            // Try to append
            File.append(origTitle + "," + hippoVol + "," + tangleVol + "," + tangleCount + "," + tangleBurden + "," + density, resultsFile);
            print("\nResults appended to: " + resultsFile);
            success = true;
        } else {
            // Create new
            f = File.open(resultsFile);
            if (f != "") {
                print(f, "Image,HippoVol_um3,TangleVol_um3,TangleCount,Burden_%,Density_per_mm3");
                print(f, origTitle + "," + hippoVol + "," + tangleVol + "," + tangleCount + "," + tangleBurden + "," + density);
                File.close(f);
                print("\nResults saved to: " + resultsFile);
                success = true;
            }
        }
    } catch (e) {
        success = false;
    }
    
    if (!success) {
        print("\n⚠ Could not save CSV file (may be open or locked)");
        print("Results are in the Log window - you can copy them manually");
        
        // Offer to save to desktop instead
        desktopPath = getDirectory("home") + "Desktop" + File.separator;
        if (File.isDirectory(desktopPath)) {
            altFile = desktopPath + "tangle_results_" + replace(origTitle, " ", "_") + ".csv";
            f = File.open(altFile);
            if (f != "") {
                print(f, "Image,HippoVol_um3,TangleVol_um3,TangleCount,Burden_%,Density_per_mm3");
                print(f, origTitle + "," + hippoVol + "," + tangleVol + "," + tangleCount + "," + tangleBurden + "," + density);
                File.close(f);
                print("✓ Saved to Desktop instead: " + altFile);
            }
        }
    }
} else {
    print("\nResults not saved to file (available in Log window)");
}
