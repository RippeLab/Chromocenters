############################################### Plot CC ######################################################################
# Creates different types of boxplots from a curated segmentation dataset generated using curateCC 
# Boxplots will be exported as PDF and stored in the folder _plotCC within the curateCC_results folder
# User can define how to filter for good/bad segmentation quality of cells for plotting

####################################################################################################################
############################################### ADJUSTABLE PARAMETERS ##############################################
####################################################################################################################

# specify folder location of the curateCC results
folder_curated <- "/Users/lfrank/Desktop/PCH-dCas9_working_copy/test_images/dCas9-EGFP/2020-06-26_16-57-28_segmentCC_results/2020-06-30_10-39-01_curateCC_results/"

# select cell quality filtering 
# "g" (include good cells only), "b" (include bad cells only), "gb" (include both)
cell_filter <- "g"

####################################################################################################################
############################################### FIXED ##############################################################

# create timestamp variable for this run
date_time <- gsub(" ","_",gsub(":","-",substr(Sys.time(),1,19)))

# create folder to export plot PDFs
dir.create(paste0(folder_curated,date_time,"_plotCC_PDFs"))

# retrieve curated area & intensity quantitifaction data from file
data_curated <- read.table(file=paste0(folder_curated,list.files(folder_curated)[grep("mean_sd_projection_curated",list.files(folder_curated))]), header=T, stringsAsFactors = F)

# check if annotation column contains unwanted characters
if(length(table(data_curated$annotation))!=3 & length(table(data_curated$annotation))<2){
  print("WARNING: The $annotation column of your dataset seems to be incomplete or contain more levels than g/b/n. This may cause errors when filtering cells based on the annotation parameter.")
  fillgaps <- readline('Set cells without g/b/n annotation to n(none) ? yes(y) or no(n):      ')
  if(fillgaps=="y"){
    gaps <- c(which(data_curated$annotation!="g" & data_curated$annotation!="b" & data_curated$annotation!="n"))
    for(n in 1:length(gaps)){
      data_curated[gaps[n],25] <- "n"
    }
  }else{
    print("NOTE:Proceeded without replacement. This may cause errors during plotting.")
  }
  }

# select cells based on user-defined filter strategy
if (cell_filter=="gb"){
  data_curated_filtered <- data_curated[data_curated[,"annotation"]=="g"|data_curated[,"annotation"]=="b",]
  if(length(rownames(data_curated_filtered))==0){stop("No cells with specified quality annotation were found. Check filter settings and whether the dataset contains cells with this type of annotation.")}
}else{
  data_curated_filtered <- data_curated[data_curated[,"annotation"]==cell_filter,]
  if(length(rownames(data_curated_filtered))==0){stop("No cells with specified quality annotation were found. Check filter settings and whether the dataset contains cells with this type of annotation.")}
}

############################### generate boxplots as PDFs #################################################
line_thickness <- 2

# boxplots for average absolute nuclear intensity values
pdf(paste0(folder_curated,date_time,"_plotCC_PDFs","/dapi_dcas9_ac_whole_nucleus.pdf"), height=5, width=4)
boxplot(values ~ vars, data = data.frame(values = c(data_curated_filtered[,"mean_dapi"],data_curated_filtered[,"mean_dcas9"],data_curated_filtered[,"mean_ac"]), vars = rep(c("a","b","c"), times=c(length(data_curated_filtered[,"mean_dapi"]),length(data_curated_filtered[,"mean_dcas9"]),length(data_curated_filtered[,"mean_ac"])))), names=c("DAPI", "dCas9", "H3K27ac"), main="Average nuclear levels", ylab="Mean nuclear intensity (a.u.)", boxlwd=line_thickness, staplelwd=line_thickness, whisklwd=line_thickness, outlwd=line_thickness, whisklty=1, staplewex=.2, boxfill="gray", las=1, outline=F)
box(lwd=line_thickness)
dev.off()

# boxplots for mean intensity enrichment in chromocenter mask versus the nucleoplasm mask
pdf(paste0(folder_curated,date_time,"_plotCC_PDFs","/dapi_dcas9_ac_enrichment_in_cc_over_np.pdf"), height=5, width=4)
boxplot(values ~ vars, data = data.frame(values = c(data_curated_filtered[,"mean_dapi_cc"]/data_curated_filtered[,"mean_dapi_np"],data_curated_filtered[,"mean_dcas9_cc"]/data_curated_filtered[,"mean_dcas9_np"],data_curated_filtered[,"mean_ac_cc"]/data_curated_filtered[,"mean_ac_np"]), vars = rep(c("a","b","c"), times=c(length(data_curated_filtered[,"mean_dapi_cc"]/data_curated_filtered[,"mean_dapi_np"]),length(data_curated_filtered[,"mean_dcas9_cc"]/data_curated_filtered[,"mean_dcas9_np"]),length(data_curated_filtered[,"mean_ac_cc"]/data_curated_filtered[,"mean_ac_np"])))), names=c("DAPI", "dCas9", "H3K27ac"), main="Enrichment in CC over NP", ylab="Normalized intensity (CC/NP)", boxlwd=line_thickness, staplelwd=line_thickness, whisklwd=line_thickness, outlwd=line_thickness, whisklty=1, staplewex=.2, boxfill="gray", las=1, outline=F)
box(lwd=line_thickness)
dev.off()

# boxplot depicting the area occupied by dCas9-segmented chromocenters relative to the total area of the nucleus
pdf(paste0(folder_curated,date_time,"_plotCC_PDFs","/dcas9_area_cc.pdf"), height=5, width=2.5)
boxplot(values ~ vars, data = data.frame(values = c(data_curated_filtered[,"area_cc"]/data_curated_filtered[,"area"]*100), vars = rep(c("a"), times=c(length(data_curated_filtered[,"area_cc"]/data_curated_filtered[,"area"])))), names=c("area"), main="Chromocenter area", ylab="Nuclear area occupied by CC (%)", boxlwd=line_thickness, staplelwd=line_thickness, whisklwd=line_thickness, outlwd=line_thickness, whisklty=1, staplewex=.2, boxfill="gray", las=1, outline=F, ylim=c(0,50))
box(lwd=line_thickness)
dev.off()

###########################################################################################################
