############################################ Segment CC ############################################################
# Thresholding based segmentation strategies to generate nucleus, chromocenter and nucleoplasm masks
# Intensities and pixel area are quantified for the masks to extract relevant image features.  

# load libraries and functions
library(EBImage)
library(abind)
source("./makeChromocenterMask.R")
source("./makeNucMask.R")

####################################################################################################################
############################################### ADJUSTABLE PARAMETERS ##############################################
####################################################################################################################

# location of the maximum-intensity .tif files (input files)
folder <- "/Users/lfrank/Desktop/PCH-dCas9_working_copy/test_images/dCas9-EGFP/"

# adjustable adaptive thresholding offset for nucleus segmentation
# default is 0.0008
threshold <- 0.0008

# select the channel in which the nuclei should be segmented
nucleus_seg_channel <- "dapi" #"ac", "dcas9", "dapi"

# remove incomplete nuclei at the image borders? "y" = yes, "n" = no
remove_bordercells <- "y"

# thresholding parameter for CC segmentation in the dCas9 channel using the makeChromocenterMask fuction
# the intensity threshold above which a pixel is considered to belong to a chromocenter (based on dCas9 enrichment)
# is given by threshold = median (nuclear intensity) +sat_cutoff*(maximum(nuclear intensity)-median(nuclear intensity))
# default is 0.1
sat_cutoff <- 0.15 

# CC segmentation will only be carried out on cells that express dCas9 (positive cells)
# default is 1.5
# A cell is considered positive if its mean nuclear dCas9 intensity is 1.5-fold above the median (background) of the total image intensity
# This parameter heavily depends on the type of construct used and the signal-to-background ratio 
mean_cutoff <- 1.5

# size-filter for nucleus segmentation (pixels), needs to be adjusted to cell size in sample
# default is 8000
nucSize_cutoff <- 8000 

# maximum number of nuclei expected per image; required for variable initialization
# default is 20
mn <- 20 

####################################################################################################################
############################################### FIXED ##############################################################

# create timestamp variable for this run
date_time <- gsub(" ","_",gsub(":","-",substr(Sys.time(),1,19)))

# create results folders
dir.create(paste0(folder,date_time,"_segmentCC_results"))
results_folder <- paste0(folder,date_time,"_segmentCC_results")
dir.create(paste0(results_folder,"/",date_time,"_nuc_shape_features"))
  
# document adjustable segmentation parameters used in this run
parameters_df<- data.frame(threshold, nucleus_seg_channel, remove_bordercells, sat_cutoff, mean_cutoff, nucSize_cutoff, mn)
colnames(parameters_df) <- c("threshold","nucleus_seg_channel","remove_bordercells","sat_cutoff","mean_cutoff","nucSize_cutoff", "mn")
write.table(file=paste0(results_folder,"/",date_time,"_segmentation_parameters.csv"), parameters_df, sep = "\t", append = F, row.names = F, quote = F)

# generate a position list to document the files that were processed in this run
positions <- list.files(path=folder, pattern="tif")

# write position list
write.table(file=paste0(results_folder,"/",date_time,"_positions.csv"), cbind(seq(1,length(positions),by=1) ,positions), sep = "\t", append = F, row.names = F, quote = F)

# initialization of variables and empty vectors for storing quantification results
mean_ac_intensity <- rep(0, mn*length(positions))
sd_ac_intensity <- rep(0, mn*length(positions))
mean_dapi_intensity <- rep(0, mn*length(positions))
sd_dapi_intensity <- rep(0, mn*length(positions))
mean_dcas9_intensity <- rep(0, mn*length(positions))
sd_dcas9_intensity <- rep(0, mn*length(positions))
area <- rep(0, mn*length(positions))
pos <- rep(0, mn*length(positions))
mean_ac_cc_intensity <- rep(0, mn*length(positions))
sd_ac_cc_intensity <- rep(0, mn*length(positions))
mean_ac_np_intensity <- rep(0, mn*length(positions))
sd_ac_np_intensity <- rep(0, mn*length(positions))
mean_dapi_cc_intensity <- rep(0, mn*length(positions))
sd_dapi_cc_intensity <- rep(0, mn*length(positions))
mean_dapi_np_intensity <- rep(0, mn*length(positions))
sd_dapi_np_intensity <- rep(0, mn*length(positions))
mean_dcas9_cc_intensity <- rep(0, mn*length(positions))
sd_dcas9_cc_intensity <- rep(0, mn*length(positions))
mean_dcas9_np_intensity <- rep(0, mn*length(positions))
sd_dcas9_np_intensity <- rep(0, mn*length(positions))
area_cc <- rep(0, mn*length(positions))
area_np <- rep(0, mn*length(positions))
num_cc <- rep(0, mn*length(positions))
cnt <- 1

# analyze nuclei
for(i in 1:length(positions)) {
  data <- readImage(paste0(folder,"/",positions[i])) #
  
  # create image objects from data
  dcas9_projection <- data[,,seq(1,numberOfFrames(data)/3,by=1)] #channel pos 1: dCas9-EGFP-VPR
  dapi_projection <- data[,,seq(numberOfFrames(data)/3+1,2*numberOfFrames(data)/3,by=1)] #channel pos 2: DAPI
  ac_projection <- data[,,seq(2*numberOfFrames(data)/3+1,3*numberOfFrames(data)/3,by=1)] #channel pos 3: H3K27ac
  
  # generate nucleus masks using the above-defined threshold offset and channel selection
  if(nucleus_seg_channel=="dapi"){nucmask <- makeNucMask(gblur(dapi_projection, sigma=3), threshold, nucSize_cutoff, remove_bordercells)}
  if(nucleus_seg_channel=="ac"){nucmask <- makeNucMask(gblur(ac_projection, sigma=3), threshold, nucSize_cutoff, remove_bordercells)}
  if(nucleus_seg_channel=="dcas9"){nucmask <- makeNucMask(gblur(dcas9_projection, sigma=3), threshold, nucSize_cutoff, remove_bordercells)}
  
  # initialize variables for CC segmentation 
  allcc <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
  allnp <- matrix(0, nrow=dim(nucmask)[1], ncol=dim(nucmask)[2])
  
  # if nuclear masks found, quantify intensities, SD and area in each mask throughout all channels
  if(length(table(nucmask))>1) {
    rois <- as.numeric(names(table(nucmask)[2:length(table(nucmask))]))
    
    for(j in 1:length(rois)) {
      mean_ac_intensity[cnt] <- mean(ac_projection[nucmask==rois[j]])
      sd_ac_intensity[cnt] <- sd(ac_projection[nucmask==rois[j]])
      mean_dapi_intensity[cnt] <- mean(dapi_projection[nucmask==rois[j]])
      sd_dapi_intensity[cnt] <- sd(dapi_projection[nucmask==rois[j]])
      mean_dcas9_intensity[cnt] <- mean(dcas9_projection[nucmask==rois[j]])
      sd_dcas9_intensity[cnt] <- sd(dcas9_projection[nucmask==rois[j]])
      area[cnt] <- sum(nucmask==rois[j])
      pos[cnt] <- i
      
      # compute nucleus shape features and positions, export them to a file (they are required for later manual curation)
      nuc_shape_features <-cbind(i,j,computeFeatures.moment(nucmask==rois[j]),computeFeatures.shape(nucmask==rois[j]))
      write.table(file=paste0(paste0(results_folder,"/",date_time,"_nuc_shape_features"),"/nuc_shape_features_pos_",i,"_cell_",j), nuc_shape_features, sep = "\t", append = F, row.names = F, quote = F)
      
      # identify if current nucleus is dCas9-positive
      # if yes, proceed with chromocenter segmentation for this nucleus
      if(mean(gblur(dcas9_projection,sigma=1)[nucmask==rois[j]])>mean_cutoff*median(gblur(dcas9_projection,sigma=1))) { # transfected cells 
        
        # generate CC masks
        mask <- nucmask
        mask[mask!=rois[j]] <- 0
        nuc <- gblur(dcas9_projection, sigma=1)
        nuc[mask!=rois[j]] <- NA
        ccmask <- makeChromocenterMask(nuc, median(nuc, na.rm=T)+sat_cutoff*(max(nuc, na.rm=T)-median(nuc, na.rm=T)), mask)
        ccmask[nuc<(median(nuc, na.rm=T)+sat_cutoff*(max(nuc, na.rm=T)-median(nuc, na.rm=T)))] <- 0
        num_cc[cnt] <- length(table(ccmask))-1
        ccmask[ccmask!=0] <- 1
        allcc[ccmask>0] <- 1
        
        # make mask for nucleoplasm in the current nucleus which is the inverse of the CC masks
        npmask <- nuc<(median(nuc, na.rm=T)+sat_cutoff*(max(nuc, na.rm=T)-median(nuc, na.rm=T)))
        npmask[is.na(npmask)] <- 0
        allnp[npmask>0] <- 1
          
        # CC region and nucleoplasm quantification
        mean_ac_cc_intensity[cnt] <- mean(ac_projection[ccmask>0])
        sd_ac_cc_intensity[cnt] <- sd(ac_projection[ccmask>0])
        mean_ac_np_intensity[cnt] <- mean(ac_projection[npmask>0])
        sd_ac_np_intensity[cnt] <- sd(ac_projection[npmask>0])
        mean_dapi_cc_intensity[cnt] <- mean(dapi_projection[ccmask>0])
        sd_dapi_cc_intensity[cnt] <- sd(dapi_projection[ccmask>0])
        mean_dapi_np_intensity[cnt] <- mean(dapi_projection[npmask>0])
        sd_dapi_np_intensity[cnt] <- sd(dapi_projection[npmask>0])
        mean_dcas9_cc_intensity[cnt] <- mean(dcas9_projection[ccmask>0])
        sd_dcas9_cc_intensity[cnt] <- sd(dcas9_projection[ccmask>0])
        mean_dcas9_np_intensity[cnt] <- mean(dcas9_projection[npmask>0])
        sd_dcas9_np_intensity[cnt] <- sd(dcas9_projection[npmask>0])
        area_cc[cnt] <- sum(ccmask>0)
        area_np[cnt] <- sum(npmask>0)
      }
      # increase total nucleus counter
      cnt <- cnt + 1
    }
  }
  
  # write out images with nucleus and CC masks overlayed onto DAPI/dCas9/ac for visualization & later curation
  # nucleus mask: yellow outlines, CC mask: blue outlines, nucleoplasm: opaque red area
  dapi_projection_with_mask <- toRGB(dapi_projection/max(dapi_projection))
  dapi_projection_with_mask <- paintObjects(nucmask, dapi_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  dapi_projection_with_mask <- paintObjects(allcc, dapi_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  dapi_projection_with_mask <- paintObjects(allnp, dapi_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(dapi_projection_with_mask, toRGB(dapi_projection/max(dapi_projection))), bits.per.sample=8, files=paste0(results_folder,"/",date_time,"_",gsub(".tif","_dapi_with_masks.tif",positions[i])))

  dcas9_projection_with_mask <- toRGB(dcas9_projection/max(dcas9_projection))
  dcas9_projection_with_mask <- paintObjects(nucmask, dcas9_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  dcas9_projection_with_mask <- paintObjects(allcc, dcas9_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  dcas9_projection_with_mask <- paintObjects(allnp, dcas9_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(dcas9_projection_with_mask, toRGB(dcas9_projection/max(dcas9_projection))), bits.per.sample=8, files=paste0(results_folder,"/",date_time,"_",gsub(".tif","_dcas9_with_masks.tif",positions[i])))
  
  ac_projection_with_mask <- toRGB(ac_projection/max(ac_projection))
  ac_projection_with_mask <- paintObjects(nucmask, ac_projection_with_mask, opac=c(1,0), col=c("yellow","yellow"))
  ac_projection_with_mask <- paintObjects(allcc, ac_projection_with_mask, opac=c(1,0), col=c("blue","blue"))
  ac_projection_with_mask <- paintObjects(allnp, ac_projection_with_mask, opac=c(0,0.2), col=c("red","red"))
  writeImage(combine(ac_projection_with_mask, toRGB(ac_projection/max(ac_projection))), bits.per.sample=8, files=paste0(results_folder,"/",date_time,"_",gsub(".tif","_ac_with_masks.tif",positions[i])))
  
  # introduce position dependent cell numbering
  cell_number_per_pos <- vector()
  pos <- pos[1:length(grep(FALSE,pos[]==0))]
  for (t in 1:max(pos)){
    cell_number_per_pos <- c(cell_number_per_pos, c(1:table(pos)[[t]]))
  }
  
}

## collect nucleus shape and position features
nuc_shape_folder <- paste0(results_folder,"/",list.files(results_folder)[grep("nuc_shape_features",list.files(results_folder))])
nuc_shape_files <- list.files(nuc_shape_folder)

nucleus_features <- data.frame()
for(b in 1:length(nuc_shape_files)){
  nucleus_features[b,1:length(colnames(read.table(file=paste0(nuc_shape_folder,"/",nuc_shape_files[b]),header=T)))] <- read.table(file=paste0(nuc_shape_folder,"/",nuc_shape_files[b]),header=T,)
}

#re-order nucleus_features according to columns position and cell (i and j) and output to results file
nucleus_features <- nucleus_features[order(nucleus_features$i,nucleus_features$j),]
colnames(nucleus_features)[1] <- "position"
colnames(nucleus_features)[2] <- "cell"
write.table(file=paste0(results_folder,"/",date_time, "_nucleus_shape_features.csv"), nucleus_features, sep = "\t", append = F, row.names = F, quote = F)
if (file.exists(nuc_shape_folder)){unlink(nuc_shape_folder, recursive = T)}

# write position/cell indices,  intensity and area results into result file
res <- cbind(pos[1:(cnt-1)],cell_number_per_pos,mean_dapi_intensity[1:(cnt-1)],sd_dapi_intensity[1:(cnt-1)],mean_ac_intensity[1:(cnt-1)],sd_ac_intensity[1:(cnt-1)],mean_dcas9_intensity[1:(cnt-1)],sd_dcas9_intensity[1:(cnt-1)],area[1:(cnt-1)],mean_dapi_cc_intensity[1:(cnt-1)],sd_dapi_cc_intensity[1:(cnt-1)],mean_dapi_np_intensity[1:(cnt-1)],sd_dapi_np_intensity[1:(cnt-1)],mean_ac_cc_intensity[1:(cnt-1)],sd_ac_cc_intensity[1:(cnt-1)],mean_ac_np_intensity[1:(cnt-1)],sd_ac_np_intensity[1:(cnt-1)],mean_dcas9_cc_intensity[1:(cnt-1)],sd_dcas9_cc_intensity[1:(cnt-1)],mean_dcas9_np_intensity[1:(cnt-1)],sd_dcas9_np_intensity[1:(cnt-1)],area_cc[1:(cnt-1)],area_np[1:(cnt-1)],num_cc[1:(cnt-1)])
colnames(res) <- c("position","cell", "mean_dapi", "sd_dapi", "mean_ac", "sd_ac", "mean_dcas9", "sd_dcas9", "area", "mean_dapi_cc", "sd_dapi_cc", "mean_dapi_np", "sd_dapi_np", "mean_ac_cc", "sd_ac_cc", "mean_ac_np", "sd_ac_np", "mean_dcas9_cc", "sd_dcas9_cc", "mean_dcas9_np", "sd_dcas9_np",  "area_cc", "area_np", "num_cc")
write.table(file=paste0(results_folder,"/",date_time, "_mean_sd_projection.csv"), res, sep = "\t", append = F, row.names = F, quote = F)


