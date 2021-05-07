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
  * exemplary boxplots as PDFs depicting chromocenter area and fold-enrichment of signals in chromocenters

## How is this script used?

### Requirements
* R version 3.6 or higher
* R packages:
  * [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html) version 4.26

### Running the analysis on the test data set
* Install the listed requirements via [Bioconductor](https://bioconductor.org/) (for [EBImage](https://bioconductor.org/packages/release/bioc/html/EBImage.html))
* Download the functions [makeNucMask](https://github.com/RippeLab/ChromocenterFeatures/blob/main/makeNucMask.R), [makeChromocenterMask](https://github.com/RippeLab/ChromocenterFeatures/blob/main/makeChromocenterMask.R), the [example data](https://github.com/RippeLab/ChromocenterFeatures/blob/main/sample_data.zip) as well as the R scripts [SegmentCC.R](https://github.com/RippeLab/ChromocenterFeatures/blob/main/segmentCC.R), [CurateCC.R](https://github.com/RippeLab/ChromocenterFeatures/blob/main/curateCC.R) and [PlotCC.R](https://github.com/RippeLab/ChromocenterFeatures/blob/main/plotCC.R) into one folder on your local computer.
* Run the R script SegmentCC in R Studio 
* Open the R script CurateCC in R Studio, define the location of the SegmentCC results folder and run the script. Follow the instructions for inspecting and scoring the segmented images.
* Open the R script PlotCC in R Studio, define the location of the CurateCC results folder and run the script to generate exemplary plots. 

## Associated scientific publications
This is associated with the following publications:
* Frank L, Weinmann R, Erdel F, Trojanowski J and Rippe K (2021) Transcriptional Activation of Heterochromatin by Recruitment of dCas9 Activators. Enhancers and Promoters - Methods and Protocols. *In Press*.
* Erdel F, Rademacher A, Vlijm R, Tünnermann J, Frank L, Weinmann R, Schweigert E, Yserentant K, Hummert J, Bauer C, Schumacher S, Al Alwash A, Normand C, Herten DP, Engelhardt J, Rippe K. Mouse Heterochromatin Adopts Digital Compaction States without Showing Hallmarks of HP1-Driven Liquid-Liquid Phase Separation. *Mol Cell*. 2020 Apr 16;78(2):236-249.e7. PMID: 3210700. PMCID: [PMC7163299](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc7163299/). DOI: [10.1016/j.molcel.2020.02.005](https://doi.org/10.1016/j.molcel.2020.02.005).
