// ===== MULTI-SAMPLE TANGLE ANALYSIS STATISTICS =====
// Analyzes the Tangle_Analysis_Summary table across all samples
// Generates group statistics, comparisons, and visualizations

print("\\Clear");
print("=================================");
print("MULTI-SAMPLE ANALYSIS");
print("=================================");

// Check if summary table exists
if (!isOpen("Tangle_Analysis_Summary")) {
    print("ERROR: Tangle_Analysis_Summary table not found!");
    print("Please run the post-processing script on your samples first.");
    exit();
}

selectWindow("Tangle_Analysis_Summary");
n = Table.size;

if (n == 0) {
    print("ERROR: Summary table is empty!");
    exit();
}

print("Analyzing " + n + " samples...\n");

// ====================================================================
// READ DATA FROM SUMMARY TABLE
// ====================================================================
samples = newArray(n);
hippoVols = newArray(n);
tangleCounts = newArray(n);
burdens = newArray(n);
densities = newArray(n);
medianVols = newArray(n);
meanVols = newArray(n);
sphericities = newArray(n);

for (i = 0; i < n; i++) {
    samples[i] = Table.getString("Sample", i);
    hippoVols[i] = Table.get("Hippo_Vol_um3", i);
    tangleCounts[i] = Table.get("Tangle_Count", i);
    burdens[i] = Table.get("Burden_%", i);
    densities[i] = Table.get("Density_per_mm3", i);
    medianVols[i] = Table.get("Median_Vol_um3", i);
    meanVols[i] = Table.get("Mean_Vol_um3", i);
    sphericities[i] = Table.get("Mean_Sphericity", i);
}

// ====================================================================
// CALCULATE GROUP STATISTICS
// ====================================================================
print("=== HIPPOCAMPUS VOLUME ===");
Array.getStatistics(hippoVols, minH, maxH, meanH, stdH);
sortedH = Array.copy(hippoVols);
Array.sort(sortedH);
medianH = sortedH[floor(n/2)];

print("Mean: " + d2s(meanH, 0) + " ± " + d2s(stdH, 0) + " µm³");
print("Median: " + d2s(medianH, 0) + " µm³");
print("Range: " + d2s(minH, 0) + " - " + d2s(maxH, 0) + " µm³");
print("CV: " + d2s(100*stdH/meanH, 1) + " %");

print("\n=== TANGLE COUNT ===");
Array.getStatistics(tangleCounts, minTC, maxTC, meanTC, stdTC);
sortedTC = Array.copy(tangleCounts);
Array.sort(sortedTC);
medianTC = sortedTC[floor(n/2)];

print("Mean: " + d2s(meanTC, 1) + " ± " + d2s(stdTC, 1) + " tangles");
print("Median: " + d2s(medianTC, 0) + " tangles");
print("Range: " + d2s(minTC, 0) + " - " + d2s(maxTC, 0));
print("CV: " + d2s(100*stdTC/meanTC, 1) + " %");

print("\n=== TANGLE BURDEN ===");
Array.getStatistics(burdens, minB, maxB, meanB, stdB);
sortedB = Array.copy(burdens);
Array.sort(sortedB);
medianB = sortedB[floor(n/2)];

print("Mean: " + d2s(meanB, 4) + " ± " + d2s(stdB, 4) + " %");
print("Median: " + d2s(medianB, 4) + " %");
print("Range: " + d2s(minB, 4) + " - " + d2s(maxB, 4) + " %");
print("CV: " + d2s(100*stdB/meanB, 1) + " %");

print("\n=== TANGLE DENSITY ===");
Array.getStatistics(densities, minD, maxD, meanD, stdD);
sortedD = Array.copy(densities);
Array.sort(sortedD);
medianD = sortedD[floor(n/2)];

print("Mean: " + d2s(meanD, 0) + " ± " + d2s(stdD, 0) + " tangles/mm³");
print("Median: " + d2s(medianD, 0) + " tangles/mm³");
print("Range: " + d2s(minD, 0) + " - " + d2s(maxD, 0));
print("CV: " + d2s(100*stdD/meanD, 1) + " %");

print("\n=== TANGLE SIZE ===");
Array.getStatistics(meanVols, minMV, maxMV, meanMV, stdMV);
sortedMV = Array.copy(meanVols);
Array.sort(sortedMV);
medianMV = sortedMV[floor(n/2)];

print("Mean across samples: " + d2s(meanMV, 2) + " ± " + d2s(stdMV, 2) + " µm³");
print("Median: " + d2s(medianMV, 2) + " µm³");
print("Range: " + d2s(minMV, 2) + " - " + d2s(maxMV, 2) + " µm³");

print("\n=== TANGLE SHAPE (Sphericity) ===");
Array.getStatistics(sphericities, minS, maxS, meanS, stdS);
print("Mean: " + d2s(meanS, 3) + " ± " + d2s(stdS, 3));
print("Range: " + d2s(minS, 3) + " - " + d2s(maxS, 3));

// ====================================================================
// CREATE STATISTICS SUMMARY TABLE
// ====================================================================
print("\n=== CREATING STATISTICS TABLE ===");

Table.create("Group_Statistics");

// Add metrics as rows
metrics = newArray("Hippocampus_Volume_um3", "Tangle_Count", "Burden_%", 
                   "Density_per_mm3", "Mean_Tangle_Size_um3", "Mean_Sphericity");
means = newArray(meanH, meanTC, meanB, meanD, meanMV, meanS);
stds = newArray(stdH, stdTC, stdB, stdD, stdMV, stdS);
medians = newArray(medianH, medianTC, medianB, medianD, medianMV, 0);
mins = newArray(minH, minTC, minB, minD, minMV, minS);
maxs = newArray(maxH, maxTC, maxB, maxD, maxMV, maxS);

for (i = 0; i < metrics.length; i++) {
    Table.set("Metric", i, metrics[i]);
    Table.set("N", i, n);
    Table.set("Mean", i, means[i]);
    Table.set("Std_Dev", i, stds[i]);
    Table.set("SEM", i, stds[i] / sqrt(n));
    Table.set("Median", i, medians[i]);
    Table.set("Min", i, mins[i]);
    Table.set("Max", i, maxs[i]);
    Table.set("CV_%", i, 100 * stds[i] / means[i]);
}

Table.update;
print("✓ Group_Statistics table created");

// ====================================================================
// VISUALIZATIONS
// ====================================================================
print("\n=== CREATING PLOTS ===");

// 1. BAR PLOT - Burden by Sample
Plot.create("Tangle Burden by Sample", "Sample", "Burden (%)");
Plot.setFontSize(12);

xPos = newArray(n);
for (i = 0; i < n; i++) {
    xPos[i] = i + 1;
}

Plot.setColor("blue");
Plot.add("bar", xPos, burdens);

// Add mean line
Plot.setColor("red");
Plot.setLineWidth(2);
Plot.drawLine(0, meanB, n + 1, meanB);
Plot.addText("Mean: " + d2s(meanB, 4) + "%", n * 0.7, meanB * 1.1);

Plot.setLimits(0, n + 1, 0, maxB * 1.2);
Plot.show();

print("✓ Burden bar plot created");

// 2. BAR PLOT - Density by Sample
Plot.create("Tangle Density by Sample", "Sample", "Density (tangles/mm³)");
Plot.setFontSize(12);

Plot.setColor("green");
Plot.add("bar", xPos, densities);

Plot.setColor("red");
Plot.setLineWidth(2);
Plot.drawLine(0, meanD, n + 1, meanD);
Plot.addText("Mean: " + d2s(meanD, 0), n * 0.7, meanD * 1.1);

Plot.setLimits(0, n + 1, 0, maxD * 1.2);
Plot.show();

print("✓ Density bar plot created");

// 3. SCATTER - Burden vs Density
Plot.create("Burden vs Density", "Burden (%)", "Density (tangles/mm³)");
Plot.setFontSize(12);

Plot.setColor("blue");
Plot.add("circle", burdens, densities);

// Add correlation
sumX = 0; sumY = 0; sumXY = 0; sumX2 = 0; sumY2 = 0;
for (i = 0; i < n; i++) {
    sumX += burdens[i];
    sumY += densities[i];
    sumXY += burdens[i] * densities[i];
    sumX2 += burdens[i] * burdens[i];
    sumY2 += densities[i] * densities[i];
}
r = (n * sumXY - sumX * sumY) / sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));

Plot.addText("r = " + d2s(r, 3), maxB * 0.7, maxD * 0.9);

Plot.setLimits(0, maxB * 1.1, 0, maxD * 1.1);
Plot.show();

print("✓ Correlation plot created");
print("  Pearson r = " + d2s(r, 3));

// 4. BOX PLOT - Burden Distribution
Plot.create("Burden Distribution Across Samples", "", "Burden (%)");
Plot.setFontSize(14);

p25B = sortedB[floor(n * 0.25)];
p75B = sortedB[floor(n * 0.75)];
p10B = sortedB[floor(n * 0.10)];
p90B = sortedB[floor(n * 0.90)];

boxLeft = 0.8;
boxRight = 1.2;
boxCenter = 1.0;

// Whiskers
Plot.setColor("black");
Plot.setLineWidth(2);
Plot.drawLine(boxCenter, p10B, boxCenter, p25B);
Plot.drawLine(boxCenter, p75B, boxCenter, p90B);
Plot.drawLine(boxLeft - 0.1, p10B, boxRight + 0.1, p10B);
Plot.drawLine(boxLeft - 0.1, p90B, boxRight + 0.1, p90B);

// Box
Plot.setColor("cyan");
xBox = newArray(boxLeft, boxRight, boxRight, boxLeft, boxLeft);
yBox = newArray(p25B, p25B, p75B, p75B, p25B);
Plot.add("filled", xBox, yBox);

Plot.setColor("black");
Plot.setLineWidth(2);
Plot.add("line", xBox, yBox);

// Median
Plot.setColor("red");
Plot.setLineWidth(4);
Plot.drawLine(boxLeft, medianB, boxRight, medianB);

// Individual points
Plot.setColor("blue");
xPoints = newArray(n);
for (i = 0; i < n; i++) {
    xPoints[i] = boxCenter + (random() - 0.5) * 0.2;
}
Plot.add("circle", xPoints, burdens);

// Labels
Plot.setColor("black");
Plot.addText("Median: " + d2s(medianB, 4) + "%", 1.3, medianB);
Plot.addText("Mean: " + d2s(meanB, 4) + "%", 1.3, meanB);
Plot.addText("n = " + n, 1.3, maxB * 0.9);

Plot.setLimits(0.5, 2, 0, maxB * 1.1);
Plot.show();

print("✓ Burden box plot created");

// 5. SCATTER - Tangle Count vs Hippocampus Volume
Plot.create("Tangle Count vs Hippocampus Volume", "Hippocampus Volume (µm³)", "Tangle Count");
Plot.setFontSize(12);

Plot.setColor("blue");
Plot.add("circle", hippoVols, tangleCounts);

// Correlation
sumX = 0; sumY = 0; sumXY = 0; sumX2 = 0; sumY2 = 0;
for (i = 0; i < n; i++) {
    sumX += hippoVols[i];
    sumY += tangleCounts[i];
    sumXY += hippoVols[i] * tangleCounts[i];
    sumX2 += hippoVols[i] * hippoVols[i];
    sumY2 += tangleCounts[i] * tangleCounts[i];
}
rCount = (n * sumXY - sumX * sumY) / sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));

Plot.addText("r = " + d2s(rCount, 3), maxH * 0.7, maxTC * 0.9);

Plot.setLimits(0, maxH * 1.1, 0, maxTC * 1.1);
Plot.show();

print("✓ Count vs Volume plot created");
print("  Pearson r = " + d2s(rCount, 3));

// ====================================================================
// OUTLIER DETECTION
// ====================================================================
print("\n=== OUTLIER DETECTION ===");
print("(Using IQR method: outside 1.5 × IQR from Q1/Q3)");

IQR_B = p75B - p25B;
lowerBound = p25B - 1.5 * IQR_B;
upperBound = p75B + 1.5 * IQR_B;

print("\nBurden outlier bounds:");
print("  Lower: " + d2s(lowerBound, 4) + " %");
print("  Upper: " + d2s(upperBound, 4) + " %");

outlierCount = 0;
for (i = 0; i < n; i++) {
    if (burdens[i] < lowerBound || burdens[i] > upperBound) {
        print("  Outlier: " + samples[i] + " = " + d2s(burdens[i], 4) + " %");
        outlierCount++;
    }
}

if (outlierCount == 0) {
    print("  No outliers detected");
}

// ====================================================================
// NORMALITY TEST (Shapiro-Wilk approximation)
// ====================================================================
print("\n=== NORMALITY CHECK ===");
print("Visual inspection recommended:");
print("  • Burden CV = " + d2s(100*stdB/meanB, 1) + " %");
if (100*stdB/meanB < 30) {
    print("    → Low variability, likely normal");
} else if (100*stdB/meanB < 50) {
    print("    → Moderate variability");
} else {
    print("    → High variability, check for outliers");
}

// ====================================================================
// SAMPLE SIZE POWER ANALYSIS
// ====================================================================
print("\n=== SAMPLE SIZE NOTES ===");
print("Current n = " + n);
print("Effect size (Cohen's d) detectable with 80% power:");

if (n >= 5) {
    // Rough estimate: d ≈ 2.8 / sqrt(n)
    detectable_d = 2.8 / sqrt(n);
    print("  d ≈ " + d2s(detectable_d, 2));
    
    if (detectable_d < 1.2) {
        print("  → Good power for large effects");
    } else if (detectable_d < 2.0) {
        print("  → Can detect medium-large effects");
    } else {
        print("  → Limited to very large effects");
        print("  → Consider n ≥ 8 for better power");
    }
} else {
    print("  → n < 5: Limited statistical power");
    print("  → Recommend n ≥ 6-8 per group");
}

// ====================================================================
// FINAL SUMMARY
// ====================================================================
print("\n=================================");
print("✓ MULTI-SAMPLE ANALYSIS COMPLETE");
print("=================================");
print("\nKey Results:");
print("  Samples analyzed: " + n);
print("  Mean burden: " + d2s(meanB, 4) + " ± " + d2s(stdB, 4) + " %");
print("  Mean density: " + d2s(meanD, 0) + " ± " + d2s(stdD, 0) + " /mm³");
print("  Mean tangle size: " + d2s(meanMV, 2) + " ± " + d2s(stdMV, 2) + " µm³");
print("\nTables created:");
print("  • Group_Statistics (summary stats)");
print("  • Tangle_Analysis_Summary (raw data)");
print("\nPlots created:");
print("  • Burden by sample (bar)");
print("  • Density by sample (bar)");
print("  • Burden vs Density (scatter)");
print("  • Burden distribution (box + points)");
print("  • Count vs Volume (scatter)");
print("\nNext steps:");
print("  • File → Save As... to export tables");
print("  • Right-click plots → Save As... for figures");
print("  • Import to GraphPad/R for formal statistics");
print("=================================");