var tauWin = getString("Tau channel window title:", "C1-2xTau_S2_PIHoechst_MaxProjection.png");
var piWin  = getString("PI channel window title:",  "C3-2xTau_S2_PIHoechst_MaxProjection.png");

print("Tau window: " + tauWin);
print("PI window:  " + piWin);

run("Set Scale...", "distance=1 known=0.433 pixel=1 unit=um");
run("Set Measurements...", "area centroid count redirect=None decimal=3");

selectWindow(tauWin);
run("Duplicate...", "title=tau_analysis");
selectWindow("tau_analysis");
run("Grays");
run("8-bit");
setThreshold(200, 255);
run("Convert to Mask");
run("Analyze Particles...", "size=10-5000 circularity=0.00-1.00 show=Nothing display clear");
var nTau = nResults;
print("Tau puncta found: " + nTau);
if (nTau == 0) {
    selectWindow("tau_analysis");
    close();
    exit("No Tau puncta found.");
}
var tauX = newArray(nTau);
var tauY = newArray(nTau);
var tt   = 0;
for (tt = 0; tt < nTau; tt++) {
    tauX[tt] = getResult("X", tt);
    tauY[tt] = getResult("Y", tt);
}
run("Clear Results");
selectWindow("tau_analysis");
close();

selectWindow(piWin);
run("Duplicate...", "title=pi_analysis");
selectWindow("pi_analysis");
run("Grays");
run("8-bit");
run("Subtract Background...", "rolling=50");
run("Gaussian Blur...", "sigma=1");
setThreshold(15, 255);
run("Convert to Mask");
run("Analyze Particles...", "size=10-500 circularity=0.20-1.00 show=Nothing display clear");
var nPI = nResults;
print("PI nuclei found: " + nPI);
if (nPI == 0) {
    selectWindow("pi_analysis");
    close();
    exit("No PI nuclei found.");
}
var piX = newArray(nPI);
var piY = newArray(nPI);
var pp  = 0;
for (pp = 0; pp < nPI; pp++) {
    piX[pp] = getResult("X", pp);
    piY[pp] = getResult("Y", pp);
}
run("Clear Results");
selectWindow("pi_analysis");
close();

print("Calculating distances...");
print("PI_nucleus,Nearest_Tau_dist_um");
var pi = 0;
for (pi = 0; pi < nPI; pi++) {
    var minDist = 999999;
    var tau     = 0;
    for (tau = 0; tau < nTau; tau++) {
        var dx   = piX[pi] - tauX[tau];
        var dy   = piY[pi] - tauY[tau];
        var dist = sqrt(dx*dx + dy*dy);
        if (dist < minDist) {
            minDist = dist;
        }
    }
    print((pi+1) + "," + d2s(minDist, 2));
}
print("=== Complete ===");
