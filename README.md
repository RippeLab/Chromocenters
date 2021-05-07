# ChromocenterFeatures (SegmentCC, CurateCC, PlotCC)
Set of scripts to extract chromocenter features from fluorescence microscopy images of mouse cells.<br/>
The scripts have two associated external functions makeNucMask() and makeChromocenterMask() required for segmentation.<br/>

## What does this set of scripts do?
# SegmentCC
* segmentation of the nucleus
* segmentation of nuclear regions occupied by chromocenters (i.e. labeled with bound dCas9 or other markers)
* quantification of:
  * nucleus shape features and xy positions of nuclei within each image
  * chromocenter number & area
  * fluorescent signal intensities in chromocenters and the surrounding nucleoplasm
# CurateCC
* uses results from SegmentCC as input
* semi-automated workflow for visual inspection and curation of the data after segmentation
# PlotCC
* uses results from CurateCC as input
* exemplary routines to visualize and plot the features extracted from the images
# Outputs:
  * text file listing the positions (images) processed in each run 
  * text files listing the above-mentioned parameters for each cell
  * images of the created nuclear and chromocenter masks for visual inspection
  * fit parameters (only in R command line)
  * optional: an image with the created nuclear and aggregate masks for visual inspection (for each cell)

## How is this script used?

### Requirements
* R version 3.6 or higher
* R packages:
  * [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html) version 4.26


### Running the analysis on the test data set
* Install the listed requirements via [Bioconductor](https://bioconductor.org/) (for [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html))
* Download the [functions](https://github.com/AnneRademacher/Optodroplets/tree/main/functions) and [example_data](https://github.com/AnneRademacher/Optodroplets/tree/main/example_data) folders as well as the R script [optodroplets_induction.R](https://github.com/AnneRademacher/Optodroplets/blob/main/optodroplets_induction.R) into one folder on your local computer (note that [example_data](https://github.com/AnneRademacher/Optodroplets/tree/main/example_data) is about 250 MB and contains the full data set from Fig. 4 of Rademacher *et al.* (2021) listed below)
* Run the R script [optodroplets_induction.R](https://github.com/AnneRademacher/Optodroplets/blob/main/optodroplets_induction.R) either in R Studio or from the R console using `source("optodroplets_induction.R")`

## Associated scientific publications
This is associated with the following publications:
* Frank L, Weinmann R, Erdel F, Trojanowski J and Rippe K (2021) Transcriptional Activation of Heterochromatin by Recruitment of dCas9 Activators. Enhancers and Promoters - Methods and Protocols. *In Press*.
* Erdel F, Rademacher A, Vlijm R, Tünnermann J, Frank L, Weinmann R, Schweigert E, Yserentant K, Hummert J, Bauer C, Schumacher S, Al Alwash A, Normand C, Herten DP, Engelhardt J, Rippe K. Mouse Heterochromatin Adopts Digital Compaction States without Showing Hallmarks of HP1-Driven Liquid-Liquid Phase Separation. *Mol Cell*. 2020 Apr 16;78(2):236-249.e7. PMID: 3210700. PMCID: [PMC7163299](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc7163299/). DOI: [10.1016/j.molcel.2020.02.005](https://doi.org/10.1016/j.molcel.2020.02.005).
