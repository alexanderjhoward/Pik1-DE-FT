// ---------------------------------------------------------------------------
// Manual (or auto-seeded then hand-corrected) cell masks for YSD microscopy
//
// Two passes over the same .lif files:
//
//   MODE = "draw"     opens each series, optionally pre-populates the ROI
//                     Manager with an automatic segmentation, and waits while
//                     you add, delete or reshape ROIs. Saves a RoiSet.zip.
//
//   MODE = "measure"  reloads those saved RoiSets and writes per_cell.csv and
//                     per_image_qc.csv in exactly the format that
//                     ysd_image_analysis.R expects.
//
// Draw once, measure as many times as you like. The RoiSets are a permanent
// record of what was counted, which is worth keeping with the manuscript.
// ---------------------------------------------------------------------------

// ---- globals -------------------------------------------------------------
// In the IJ1 macro language, assigning to a name inside a function creates a
// LOCAL variable. Accumulators written to by functions must be declared with
// `var` at the top level or the updates are silently discarded.
var cellRows = newArray(0);
var qcRows   = newArray(0);
var output   = "";
var roiDir   = "";

// ============================ CONFIG =======================================

MODE = "measure";          // "draw" first, then change to "measure"

LIF_CH_LIGAND = 1;      // ligand           (PE)
LIF_CH_EXPR   = 2;      // receptor display (AF488)
LIF_CH_BF     = 3;      // brightfield

// Which channel to display while drawing. Brightfield is the least biased
// thing to outline on, since it is independent of both measured signals.
DRAW_ON = "BF";         // "BF" or "EXP"

// Only process series whose name matches this regular expression. The .lif
// holds every ligand, but only Avr-PikA and Avr-PikC are needed, so this cuts
// the drawing work by more than half. Set to "" to process every series.
SERIES_FILTER = ".*_[ACac]([^A-Za-z].*|)$";

PRESEED       = true;   // propose automatic ROIs for you to correct
SKIP_EXISTING = true;   // don't re-open series that already have a RoiSet

// Auto-seeding parameters (ignored if PRESEED is false)
MIN_SIZE      = 5;      // um^2
MAX_SIZE      = 60;     // um^2
MIN_CIRC      = 0.65;
BLUR_SIGMA    = 1.5;
THRESH_METHOD = "Otsu";
DO_WATERSHED  = true;

// Shrink every ROI by this many pixels before measuring, applied uniformly.
// Boundary pixels are a blur of cell and background, so excluding them
// removes a bias that scales with cell brightness. At 0.18 um/px, 3 px is
// about 0.5 um. Set to 0 to measure the ROI exactly as drawn.
ROI_SHRINK_PX = 3;

RECURSE = true;

// ===========================================================================

root   = getDirectory("Choose the parent folder containing your Run folders");
roiDir = root + "ROIs" + File.separator;      // assigns the global above
output = root + "Quantification" + File.separator;
File.makeDirectory(roiDir);
File.makeDirectory(output);

cellRows = Array.concat(cellRows,
  "folder,base,cell,area,expr_mean,expr_median,expr_max,expr_clip," +
  "ligand_mean,ligand_median,ligand_max,ligand_clip,expr_bg,ligand_bg");
qcRows = Array.concat(qcRows,
  "folder,base,n_cells,bit_depth,pixel_width,unit,expr_bg,ligand_bg," +
  "expr_sat_frac,ligand_sat_frac,seg_mode");

run("Set Measurements...", "area mean standard min median redirect=None decimal=4");
if (MODE == "measure") setBatchMode(true);   // drawing needs windows visible

processFolder(root);

if (MODE == "measure") {
    if (cellRows.length <= 1) {
        setBatchMode(false);
        exit("No cells were accumulated. The CSVs would contain only headers, " +
             "so nothing was written. Check that the ROIs folder contains " +
             "RoiSet.zip files and that SERIES_FILTER matches their names.");
    }
    File.saveString(String.join(cellRows, "\n"), output + "per_cell.csv");
    File.saveString(String.join(qcRows, "\n"),   output + "per_image_qc.csv");
    setBatchMode(false);
    print("Measurement done: " + (cellRows.length - 1) + " cells across " +
          (qcRows.length - 1) + " fields. CSVs written to " + output);
} else {
    print("Drawing pass done. RoiSets in " + roiDir);
    print("Now set MODE = \"measure\" at the top and run again.");
}

// ===========================================================================

function processFolder(dir) {
    list = getFileList(dir);
    for (i = 0; i < list.length; i++) {
        p = dir + list[i];
        if (File.isDirectory(p)) {
            if (RECURSE && indexOf(list[i], "Quantification") < 0
                        && indexOf(list[i], "ROIs") < 0)
                processFolder(p);
            continue;
        }
        if (!endsWith(toLowerCase(list[i]), ".lif")) continue;
        processLif(dir, list[i]);
    }
}

function processLif(dir, fname) {
    path = dir + fname;
    ftag = folderTag(dir);

    run("Bio-Formats Macro Extensions");
    Ext.setId(path);
    Ext.getSeriesCount(nSeries);

    for (sr = 0; sr < nSeries; sr++) {
        Ext.setSeries(sr);
        Ext.getSeriesName(sName);

        if (SERIES_FILTER != "" && !matches(sName, SERIES_FILTER)) continue;

        base = safeName(sName);
        roiPath = roiDir + ftag + "__" + base + "_RoiSet.zip";

        if (MODE == "draw" && SKIP_EXISTING && File.exists(roiPath)) {
            print("SKIP (already drawn): " + ftag + " / " + base);
            continue;
        }
        if (MODE == "measure" && !File.exists(roiPath)) {
            print("NO ROIs for " + ftag + " / " + base + ", skipped.");
            continue;
        }

        run("Bio-Formats Importer",
            "open=[" + path + "] autoscale color_mode=Default " +
            "view=Hyperstack stack_order=XYCZT series_" + (sr + 1));
        img = getTitle();
        getDimensions(w, h, nc, ns, nf);
        if (nc < 2) { close("*"); continue; }
        if (ns > 1) {
            run("Z Project...", "projection=[Max Intensity]");
            close(img);
            img = getTitle();
        }

        bd = bitDepth();
        getPixelSize(unit, pw, ph);

        selectWindow(img);
        run("Duplicate...", "title=LIG duplicate channels=" + LIF_CH_LIGAND);
        selectWindow(img);
        run("Duplicate...", "title=EXP duplicate channels=" + LIF_CH_EXPR);
        haveBf = (nc >= LIF_CH_BF);
        if (haveBf) {
            selectWindow(img);
            run("Duplicate...", "title=BF duplicate channels=" + LIF_CH_BF);
        }
        close(img);

        selectWindow("LIG"); ligSat = fracSaturated();
        selectWindow("EXP"); expSat = fracSaturated();

        if (MODE == "draw") drawPass(ftag, base, roiPath, haveBf);
        else measurePass(ftag, base, roiPath, bd, pw, unit, expSat, ligSat);

        close("*");
        roiManager("reset");
    }
    Ext.close();
}

// ------------------------------- DRAW --------------------------------------

function drawPass(ftag, base, roiPath, haveBf) {
    roiManager("reset");

    if (PRESEED) seedRois();

    disp = DRAW_ON;
    if (disp == "BF" && !haveBf) disp = "EXP";
    selectWindow(disp);
    run("Enhance Contrast", "saturated=0.35");
    setTool("oval");   // ovals suit yeast; press the freehand icon for buds
    roiManager("Show All with labels");

    n0 = roiManager("count");
    waitForUser("Edit cell ROIs: " + ftag + " / " + base,
        "Showing " + n0 + " proposed ROIs on the " + disp + " channel.\n \n" +
        "  - Delete a bad ROI: click it in the ROI Manager, press Delete\n" +
        "  - Add a missed cell: outline it, then press T\n" +
        "  - Reshape: select it, drag the outline, then press Update\n" +
        "  - Exclude cells touching the image edge and obvious debris\n \n" +
        "Draw on brightfield where possible so the mask is independent of\n" +
        "both measured channels.\n \n" +
        "Click OK when this field is done. Click Cancel to stop entirely.");

    n1 = roiManager("count");
    if (n1 > 0) {
        roiManager("Deselect");
        roiManager("Save", roiPath);
        print("Saved " + n1 + " ROIs for " + ftag + " / " + base +
              " (seeded " + n0 + ")");
    } else {
        print("No ROIs kept for " + ftag + " / " + base + ", nothing saved.");
    }
}

function seedRois() {
    selectWindow("EXP");
    run("Duplicate...", "title=SEG");
    selectWindow("SEG");
    run("Gaussian Blur...", "sigma=" + BLUR_SIGMA);
    setAutoThreshold(THRESH_METHOD + " dark");
    setOption("BlackBackground", true);
    run("Convert to Mask");
    run("Fill Holes");
    if (DO_WATERSHED) run("Watershed");

    getPixelSize(u, pwx, phx);
    calibrated = (pwx > 0) && (u != "pixel") && (u != "pixels") && (u != "");
    opts = "size=" + MIN_SIZE + "-" + MAX_SIZE +
           " circularity=" + MIN_CIRC + "-1.00";
    if (!calibrated) opts = opts + " pixel";
    opts = opts + " exclude add";

    run("Analyze Particles...", opts);
    close("SEG");
}

// ------------------------------ MEASURE ------------------------------------

function measurePass(ftag, base, roiPath, bd, pw, unit, expSat, ligSat) {
    roiManager("reset");
    roiManager("Open", roiPath);
    nCells = roiManager("count");
    if (nCells == 0) return;

    // Background: everything outside every drawn ROI.
    expBg = 0; ligBg = 0;
    selectWindow("EXP");
    if (nCells > 1) {
        roiManager("Deselect");
        roiManager("Combine");
    } else {
        roiManager("select", 0);
    }
    // Grow before inverting so the background excludes the blurred halo
    // around each cell rather than counting it as background.
    if (ROI_SHRINK_PX > 0)
        run("Enlarge...", "enlarge=" + ROI_SHRINK_PX + " pixel");
    run("Make Inverse");
    if (selectionType() >= 0) {
        expBg = getValue("Median");
        run("Restore Selection");
        selectWindow("LIG");
        run("Restore Selection");
        ligBg = getValue("Median");
    }
    selectWindow("EXP"); run("Select None");
    selectWindow("LIG"); run("Select None");

    for (c = 0; c < nCells; c++) {
        selectWindow("EXP");
        roiManager("select", c);
        if (ROI_SHRINK_PX > 0)
            run("Enlarge...", "enlarge=-" + ROI_SHRINK_PX + " pixel");
        a    = getValue("Area");
        mE   = getValue("Mean");
        mdE  = getValue("Median");
        maxE = getValue("Max");
        clE  = roiClipFrac();

        selectWindow("LIG");
        roiManager("select", c);
        if (ROI_SHRINK_PX > 0)
            run("Enlarge...", "enlarge=-" + ROI_SHRINK_PX + " pixel");
        mL   = getValue("Mean");
        mdL  = getValue("Median");
        maxL = getValue("Max");
        clL  = roiClipFrac();

        cellRows = Array.concat(cellRows,
            ftag + "," + base + "," + (c + 1) + "," + a + "," +
            mE + "," + mdE + "," + maxE + "," + clE + "," +
            mL + "," + mdL + "," + maxL + "," + clL + "," + expBg + "," + ligBg);
    }

    qcRows = Array.concat(qcRows,
        ftag + "," + base + "," + nCells + "," + bd + "," + pw + "," +
        unit + "," + expBg + "," + ligBg + "," + expSat + "," + ligSat +
        ",manual_shrink" + ROI_SHRINK_PX);

    print("Measured " + nCells + " cells: " + ftag + " / " + base);
}

// ------------------------------ HELPERS ------------------------------------

function safeName(s) {
    s = replace(s, "[\\\\/:*?\"<>|,]", "_");
    s = replace(s, " ", "_");
    return s;
}

function folderTag(dir) {
    parts = split(dir, File.separator);
    if (parts.length == 0) return "root";
    return parts[parts.length - 1];
}

// Fraction of pixels AT the detector ceiling within the current selection.
// A cell with a handful of clipped pixels is still measurable; one that is
// clipped across most of its area is not.
function roiClipFrac() {
    b = bitDepth();
    if (b == 8) { getHistogram(v, cts, 256);   top = 255; }
    else        { getHistogram(v, cts, 65536); top = 65535; }
    total = 0;
    for (k = 0; k < cts.length; k++) total = total + cts[k];
    if (total == 0) return 0;
    return cts[cts.length - 1] / total;
}

function fracSaturated() {
    b = bitDepth();
    if (b == 8) { getHistogram(values, counts, 256);   top = 255;   nb = 256; }
    else        { getHistogram(values, counts, 65536); top = 65535; nb = 65536; }
    total = 0; hot = 0;
    for (k = 0; k < nb; k++) {
        total = total + counts[k];
        if (k >= top - 1) hot = hot + counts[k];
    }
    if (total == 0) return 0;
    return hot / total;
}
