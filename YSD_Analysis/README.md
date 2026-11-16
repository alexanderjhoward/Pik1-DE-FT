# Modeling directed evolution data with ESM-2 to phenotype novel genomic variants of interest

## Background
For this section, I analyzed YSD imaging and FACS data to confirm predictions about Pik-1 binding behavior against Avr-Pik variants of interest.

## Method
Start by running the custom ImageJ macro for quantifying the mean receptor and ligand signal intensities of outlined cells. To do so, open Fiji, then drag-and-drop the "ysd_manual_rois.ijm" file onto the Fiji toolbar to open it. 

There are two modes for the macro that can be changed in the editor under the "MODE" variable.

"draw": When the MODE variable is set to "draw", the macro will assist with drawing cell outlines and supports the outlining of additional cells in the image. Click "Run" in the editor and select the folder containing a .lif file you want to process. The macro will make two subdirectories (ROIs/ and Quantification/) to save data to. To remove any poorly drawn ROIs, click and delete them. To draw any new ROIs I would reccommend using the oval tool and pressing the "T" key to log the ROI. When finished outlining ROIs click "OK" on the ROI Macro Popup to save them under the ROIs/ subdirectory and move on to the next image. ROIs can be re-used by the macro in "measure" mode (our ROIs can be used to re-create our exact measurments from our figures, found in Zotero).

"measure": Once all images are outlined, the macro can be used to quantify per-cell intensities and background intensities for receptor and ligand signal. Click "Run" and select the folder containing your ROIs/ subdirectory. The calculated intensitites will  get saved to the Quantification/ directory for downstream analysis.

After all per-cell intensities are collected, the values can be processed with the **ysd_image_analysis.R** script.

FACS analysis of one YSD staining session was conducted to verify differences in receptor staining in a high throughput manner. Fluorescence compensation, cell gating, and expression-normalized ligand staining were calculated using the **fcs_normalized_binding.R** script.

Visualizations of YSD receptor expression and FACS analsysis (Fig 4c and 4d) were plotted with the **Plot_YSD_FACS.ipynb** notebook.

