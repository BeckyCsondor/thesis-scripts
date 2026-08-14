// ============================================================
//  PI / Hoechst Viability Batch Macro
//  Tissue:  Mouse hippocampal organotypic slice
//  Pixel size: 0.433 um/px
//
//  Control:  Ch1=Hoechst  Ch2=PI(yellow)
//  wtTau:    Ch1=Hoechst  Ch2=Tau(red)  Ch3=PI(yellow)
//  2xTau:    Ch1=Tau(red) Ch2=Hoechst   Ch3=PI(yellow)
// ============================================================

// Parameters - nucleus sizes now in pixels² based on test results
// PI nuclei average ~34 px² (pyknotic/shrunken dead nuclei)
// Hoechst nuclei larger ~100-500 px²
var MIN_HOECHST_SIZE = 50;
var MAX_HOECHST_SIZE = 500;
var MIN_PI_SIZE      = 10;
var MAX_PI_SIZE      = 500;
var PI_THRESHOLD_LOW = 15;   // manual threshold lower bound for PI channel (confirmed final value)
var BG_RADIUS        = 50;
var THRESHOLD_HOECHST = "Otsu";
var MIN_CIRCULARITY   = 0.20;
var MAX_CIRCULARITY   = 1.00;

if (nImages == 0)
    exit("No images open. Please open your PNG files in Fiji first.");

print("=== PI/Hoechst Batch Analysis ===");
print("Images open: " + nImages);

// Snapshot original titles before loop
var nOpen      = nImages;
var origTitles = newArray(nOpen);
var t          = 0;
for (t = 0; t < nOpen; t++) {
    selectImage(t + 1);
    origTitles[t] = getTitle();
    print("Found: " + origTitles[t]);
}

// CSV setup
var csvDir = getDirectory("Choose folder to save CSV");
if (lengthOf(csvDir) == 0) exit("No folder chosen.");
var csvFile = getString("CSV filename:", "PI_Hoechst_results.csv");
if (lengthOf(csvFile) == 0) exit("No filename.");
var csvPath = csvDir + csvFile;

var fileExists = File.exists(csvPath);
if (!fileExists) {
    var hdr = "Sample,Condition,Hoechst_Ch,PI_Ch,Total_Nuclei,Dead_Nuclei,Alive_Nuclei,Seg_Pct_Viable,Seg_Pct_Dead,Hoechst_IntDen,PI_IntDen,Total_IntDen,Int_Pct_Viable,Int_Pct_Dead,Threshold_Hoechst,Threshold_PI\n";
    File.saveString(hdr, csvPath);
    print("CSV created: " + csvPath);
} else {
    print("Appending to: " + csvPath);
}

print("Processing " + nOpen + " images...");

// ============================================================
//  Helper: close all windows whose title contains a keyword.
//  Scans ALL open images by position and closes matches.
//  Does NOT use isOpen() which fails after Convert to Mask.
// ============================================================
function closeWindowsContaining(keyword) {
    var pass = 0;
    // Two passes to catch cases where closing shifts positions
    for (pass = 0; pass < 2; pass++) {
        var w = nImages;
        while (w >= 1) {
            selectImage(w);
            if (indexOf(getTitle(), keyword) >= 0) {
                close();
            }
            w = w - 1;
        }
    }
}

// ============================================================
//  MAIN BATCH LOOP
// ============================================================
var imgIdx = 0;
for (imgIdx = 0; imgIdx < nOpen; imgIdx++) {

    var origTitle = origTitles[imgIdx];
    selectWindow(origTitle);
    print("--- Processing: " + origTitle + " ---");

    // Detect condition
    var condition   = "Unknown";
    var hoechstCh   = 1;
    var piCh        = 2;
    var nChExpected = 2;
    var skipImage   = false;

    if (indexOf(origTitle, "Control") >= 0) {
        condition   = "Control";
        hoechstCh   = 1;
        piCh        = 2;
        nChExpected = 2;
    } else if (indexOf(origTitle, "wtTau") >= 0) {
        condition   = "wtTau";
        hoechstCh   = 2;
        piCh        = 3;
        nChExpected = 3;
    } else if (indexOf(origTitle, "2xTau") >= 0) {
        condition   = "2xTau";
        hoechstCh   = 2;
        piCh        = 3;
        nChExpected = 3;
    } else {
        print("WARNING: Cannot detect condition from: " + origTitle + " - skipping.");
        skipImage = true;
    }

    if (skipImage == false) {
        var nCh = nSlices;
        if (nCh != nChExpected) {
            print("WARNING: Expected " + nChExpected + " channels, found " + nCh + " - skipping.");
            skipImage = true;
        }
    }

    if (skipImage == false) {

        print("Condition: " + condition + "  Hoechst Ch: " + hoechstCh + "  PI Ch: " + piCh);

        run("Set Scale...", "distance=1 known=0.433 pixel=1 unit=um");
        run("Set Measurements...", "area mean integrated count redirect=None decimal=3");

        // Build expected window title strings
        var hoechstIntTitle = "C" + hoechstCh + "-intensity_copy";
        var piIntTitle      = "C" + piCh      + "-intensity_copy";
        var hoechstSegTitle = "C" + hoechstCh + "-seg_copy";
        var piSegTitle      = "C" + piCh      + "-seg_copy";

        // --------------------------------------------------
        //  Intensity duplicate
        // --------------------------------------------------
        selectWindow(origTitle);
        run("Duplicate...", "title=intensity_copy duplicate");
        run("Split Channels");

        // --------------------------------------------------
        //  Segmentation duplicate
        // --------------------------------------------------
        selectWindow(origTitle);
        run("Duplicate...", "title=seg_copy duplicate");
        run("Split Channels");

        // --------------------------------------------------
        //  Measure Hoechst intensity
        // --------------------------------------------------
        selectWindow(hoechstIntTitle);
        run("Grays");
        run("8-bit");
        run("Subtract Background...", "rolling=" + BG_RADIUS);
        makeRectangle(0, 0, getWidth(), getHeight());
        run("Measure");
        var hoechstIntDen = 0;
        if (nResults > 0) {
            hoechstIntDen = getResult("RawIntDen", nResults - 1);
            print("  Hoechst IntDen: " + hoechstIntDen);
        }
        run("Clear Results");

        // --------------------------------------------------
        //  Measure PI intensity
        // --------------------------------------------------
        selectWindow(piIntTitle);
        run("Grays");
        run("8-bit");
        run("Subtract Background...", "rolling=" + BG_RADIUS);
        makeRectangle(0, 0, getWidth(), getHeight());
        run("Measure");
        var piIntDen = 0;
        if (nResults > 0) {
            piIntDen = getResult("RawIntDen", nResults - 1);
            print("  PI IntDen: " + piIntDen);
        }
        run("Clear Results");

        // Close intensity windows using title-scan helper
        closeWindowsContaining("intensity_copy");

        // Intensity viability
        var totalIntDen  = hoechstIntDen + piIntDen;
        var intPctDead   = 0;
        var intPctViable = 0;
        if (totalIntDen > 0) {
            intPctDead   = (piIntDen      / totalIntDen) * 100;
            intPctViable = (hoechstIntDen / totalIntDen) * 100;
        }

        // --------------------------------------------------
        //  Segmentation - Hoechst
        // --------------------------------------------------
        selectWindow(hoechstSegTitle);
        run("Grays");
        run("8-bit");
        run("Subtract Background...", "rolling=" + BG_RADIUS);
        run("Gaussian Blur...", "sigma=1");
        setAutoThreshold(THRESHOLD_HOECHST + " dark");
        run("Convert to Mask");
        run("Fill Holes");
        run("Watershed");
        run("Set Scale...", "distance=0 known=0 unit=pixel");
        run("Analyze Particles...",
            "size=" + MIN_HOECHST_SIZE + "-" + MAX_HOECHST_SIZE +
            " circularity=" + MIN_CIRCULARITY + "-" + MAX_CIRCULARITY +
            " show=Nothing display clear summarize");

        var hoechstCount = nResults;
        run("Clear Results");

        // --------------------------------------------------
        //  Segmentation - PI
        // --------------------------------------------------
        selectWindow(piSegTitle);
        run("Grays");
        run("8-bit");
        run("Subtract Background...", "rolling=" + BG_RADIUS);
        run("Gaussian Blur...", "sigma=1");
        setThreshold(PI_THRESHOLD_LOW, 255);
        run("Convert to Mask");
        run("Fill Holes");
        run("Watershed");
        run("Set Scale...", "distance=0 known=0 unit=pixel");
        run("Analyze Particles...",
            "size=" + MIN_PI_SIZE + "-" + MAX_PI_SIZE +
            " circularity=" + MIN_CIRCULARITY + "-" + MAX_CIRCULARITY +
            " show=Nothing display clear summarize");

        var piCount = nResults;
        run("Clear Results");

        // Close seg windows using title-scan helper
        closeWindowsContaining("seg_copy");

        // Segmentation viability
        var deadCount  = piCount;
        var aliveCount = hoechstCount - deadCount;
        if (aliveCount < 0) { aliveCount = 0; }
        var segPctViable = 0;
        var segPctDead   = 0;
        if (hoechstCount > 0) {
            segPctViable = (aliveCount / hoechstCount) * 100;
            segPctDead   = (deadCount  / hoechstCount) * 100;
        }

        print("  Total nuclei (Hoechst): " + hoechstCount);
        print("  Dead nuclei  (PI+):     " + deadCount);
        print("  Seg  % Viable: " + d2s(segPctViable, 1) + " %");
        print("  Seg  % Dead:   " + d2s(segPctDead,   1) + " %");
        print("  Int  % Viable: " + d2s(intPctViable, 1) + " %");
        print("  Int  % Dead:   " + d2s(intPctDead,   1) + " %");

        // Write CSV row
        var sampleName = replace(origTitle, ".png", "");
        sampleName     = replace(sampleName, ".PNG", "");

        var r1  = sampleName + ",";
        var r2  = condition + ",";
        var r3  = hoechstCh + ",";
        var r4  = piCh + ",";
        var r5  = hoechstCount + ",";
        var r6  = deadCount + ",";
        var r7  = aliveCount + ",";
        var r8  = d2s(segPctViable, 2) + ",";
        var r9  = d2s(segPctDead,   2) + ",";
        var r10 = d2s(hoechstIntDen, 0) + ",";
        var r11 = d2s(piIntDen,      0) + ",";
        var r12 = d2s(totalIntDen,   0) + ",";
        var r13 = d2s(intPctViable,  2) + ",";
        var r14 = d2s(intPctDead,    2) + ",";
        var r15 = THRESHOLD_HOECHST + ",";
        var r16 = PI_THRESHOLD_LOW + "\n";
        var csvRow = r1+r2+r3+r4+r5+r6+r7+r8+r9+r10+r11+r12+r13+r14+r15+r16;
        File.append(csvRow, csvPath);

        selectWindow(origTitle);
        print("Done: " + origTitle);

    } // end skipImage == false

} // end loop

print("=== Batch complete. Results saved to: " + csvPath + " ===");
showMessage("Batch Complete", "All images processed.\nResults saved to:\n" + csvPath);
