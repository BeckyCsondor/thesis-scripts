// ============================================================
// mito_3D_render.ijm
// Fiji macro for 3D volume rendering of mitochondrion + granules
// ============================================================

// ---- CONFIGURATION ----
FILLED_PATH   = "/ceph/users/loo89671/mitogranule_segmentationanalysis/Pos018_mito_filled_v3.tif";
SEG_PATH      = "/ceph/users/loo89671/mitogranule_segmentationanalysis/Pos018_ImageForProc_1.2_GranulesandMembranes_Pre-trained2DU-Net_Depth5.tiff";
GRANULE_LABEL = 2;
GRAN_BLUR     = 4.0;
GRAN_THRESH   = 50;
// -----------------------

print("\\Clear");
print("=== Mito 3D render ===");

// ---- Step 1: Load filled mito volume ----
print("Loading filled mito volume...");
open(FILLED_PATH);
rename("mito_filled");
run("8-bit");
print("Mito volume loaded. Stack size: " + nSlices + " slices.");

// ---- Step 2: Load segmentation and isolate granules ----
print("Loading segmentation...");
open(SEG_PATH);
rename("segmentation");
print("Segmentation loaded. Stack size: " + nSlices + " slices.");

// Duplicate and create granule mask using thresholding
print("Extracting granule label (" + GRANULE_LABEL + ")...");
selectWindow("segmentation");
run("Duplicate...", "title=granule_mask duplicate");
selectWindow("granule_mask");

// Use setThreshold to isolate exactly the granule label value
setThreshold(GRANULE_LABEL, GRANULE_LABEL);
run("Convert to Mask", "method=Default background=Dark black");
run("8-bit");
print("Granule mask created.");

// ---- Step 3: Blur granules to smooth masses ----
print("Blurring granules (sigma=" + GRAN_BLUR + ")...");
selectWindow("granule_mask");
run("Gaussian Blur 3D...", "x=" + GRAN_BLUR + " y=" + GRAN_BLUR + " z=" + GRAN_BLUR);

// Re-threshold after blur
setThreshold(GRAN_THRESH, 255);
run("Convert to Mask", "method=Default background=Dark black");
run("8-bit");
print("Granule blur done.");

// ---- Step 4: Close raw segmentation ----
selectWindow("segmentation");
close();

// ---- Step 5: Apply LUTs for colour before 3D Viewer ----
// Set mito to blue
selectWindow("mito_filled");
run("Blue");

// Set granules to magenta
selectWindow("granule_mask");
run("Magenta");

// ---- Step 6: Open 3D Viewer ----
print("Opening 3D Viewer...");
run("3D Viewer");

// Add mito — Volume mode, translucent
selectWindow("mito_filled");
call("ij3d.ImageJ3DViewer.add", "mito_filled", "Volume", "mito_volume", "30", "true", "true", "true", "2", "0");
call("ij3d.ImageJ3DViewer.select", "mito_volume");
call("ij3d.ImageJ3DViewer.setColor", "100", "149", "237");
call("ij3d.ImageJ3DViewer.setTransparency", "0.75");

// Add granules — Volume mode, opaque
selectWindow("granule_mask");
call("ij3d.ImageJ3DViewer.add", "granule_mask", "Volume", "granules", "30", "true", "true", "true", "2", "0");
call("ij3d.ImageJ3DViewer.select", "granules");
call("ij3d.ImageJ3DViewer.setColor", "217", "70", "168");
call("ij3d.ImageJ3DViewer.setTransparency", "0.0");

print("Done! Rotate with mouse in the 3D Viewer.");
print("File > Take Snapshot to save an image.");
