###-------------------------------------------------------------------------------------------------
#####                           CREATE DATASETS FOR USE WITH FASTQTL                           #####
#
# Version: Create_mQTL_data.v1.20160404
#
# Last update: 2016-04-04
# Written by: Sander W. van der Laan (s.w.vanderlaan-2@umcutrecht.nl)
#                                                    
# Description: Script to create the necessary datasets for use with fastQTL. As per fastQTL 
#              instructions the file should look like the example below in the end. See also: 
#              http://fastqtl.sourceforge.net.
#
#              #Chr	start	end	ID	UNR1	UNR2	UNR3	UNR4 
#              chr1	173863	173864	ENSG123	-0.50 0.82	-0.71	0.83
#              chr1	685395	685396	ENSG456	-1.13	1.18	-0.03	0.11
#              chr1	700304	700305	ENSG789	-1.18	1.32	-0.36	1.26
#
#
###-------------------------------------------------------------------------------------------------
### CLEAR THE BOARD
rm(list = ls())

###-------------------------------------------------------------------------------------------------
### GENERAL R SETUP 

### FUNCTION TO INSTALL PACKAGES, VERSION A -- This is a function found by 
### Sander W. van der Laan online from @Samir: 
### http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
### Compared to VERSION 1 the advantage is that it will automatically check in both CRAN and Bioconductor
install.packages.auto <- function(x) { 
     x <- as.character(substitute(x)) 
     if (isTRUE(x %in% .packages(all.available = TRUE))) { 
          eval(parse(text = sprintf("require(\"%s\")", x)))
     } else { 
          # Update installed packages - this may mean a full upgrade of R, which in turn
          # may not be warrented. 
          #update.packages(ask = FALSE) 
          eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
     }
     if (isTRUE(x %in% .packages(all.available = TRUE))) { 
          eval(parse(text = sprintf("require(\"%s\")", x)))
     } else {
          source("http://bioconductor.org/biocLite.R")
          # Update installed packages - this may mean a full upgrade of R, which in turn
          # may not be warrented.
          #biocLite(character(), ask = FALSE) 
          eval(parse(text = sprintf("biocLite(\"%s\")", x)))
          eval(parse(text = sprintf("require(\"%s\")", x)))
     }
}
install.packages.auto("xlsx")
install.packages.auto("data.table")
install.packages.auto("dplyr")
install.packages.auto("plyr")
install.packages.auto("stringr")
install.packages.auto("stringr")

# Create datestamp
Today = format(as.Date(as.POSIXlt(Sys.time())), "%Y%m%d")

###-------------------------------------------------------------------------------------------------
### SETUP ANALYSIS
# Assess where we are
getwd()
# Set locations
METHYLDATA = "/Users/sanderwvanderlaan/PLINK/_AE_ORIGINALS/AEMethylS_IlluminaMethylation450K"
DAT_loc = paste0(METHYLDATA,"/AEM_COMPILE_DATASET/20150112_AEM_postQC7.b37.RData")

# Location of general output directory, and create it if doesn't exist
ROOT_loc = "/Users/sanderwvanderlaan/PLINK/analyses/epigenetics"
dir.create(paste0(ROOT_loc,"/mqtl"))
OUT_loc = paste0(ROOT_loc,"/mqtl")

setwd(ROOT_loc)

# Some other settings
Subset = "All" # All, Auto or X Chromosomes

### ------------------------------------------------------------------------------------------------
### LOAD AEM DATASET
cat("Loading AEMethylS covariates...")
load(DAT_loc)

### ------------------------------------------------------------------------------------------------
### PREPARE DATASET
# Let's put the originals aside
df.AEMP.original = df.AEMP
df.AEMB.original = df.AEMB

# Create a rowname variable
df.AEMP$Methylation450K_ID = rownames(df.AEMP)
df.AEMB$Methylation450K_ID = rownames(df.AEMB)
# Reset the `rownames` of your original data
rownames(df.AEMP) = NULL
rownames(df.AEMB) = NULL

# Let's edit the SampleIDs and strip "spaces"
df.AEMP$SampleID <- as.character(df.AEMP$SampleID)
df.AEMP$SampleID <- gsub("\\s" , "" , df.AEMP$SampleID)
row.names(df.AEMP) <- df.AEMP$SampleID

df.AEMB$SampleID <- as.character(df.AEMB$SampleID)
df.AEMB$SampleID <- gsub("\\s" , "" , df.AEMB$SampleID)
row.names(df.AEMB) <- df.AEMB$SampleID

###-------------------------------------------------------------------------------------------------
### PROBE FILTERING 
cat("Probe filtering...\n")
# Get a list of probes 
if (Subset == "All") {
     Probe_subset = mtx.HD450K[,1]
} else if (Subset == "Auto") {
     Probe_subset = mtx.HD450K[(mtx.HD450K$CHR != "X" & mtx.HD450K$CHR != "Y"), 1]  
} else if (Subset == "X") {
     Probe_subset = mtx.HD450K[(mtx.HD450K$CHR == "X"), 1]
}

# Get position of all CpG probes by name.
# NOTE: make sure that there are no variables in the AEDB starting with
#       - cg[0-9]
#       - ch.[0-9]
#       - ch.X
CpG_p = grep("^cg[0-9]|ch.[0-9]|ch.[X]", colnames(df.AEMP), perl = TRUE)
# clinical data for plaque methylation data
clinic_p = as.integer(labels(colnames(df.AEMP))[-CpG_p])
CpG_b = grep("^cg[0-9]|ch.[0-9]|ch.[X]", colnames(df.AEMB), perl = TRUE)
# clinical data for blood methylation data
clinic_b = as.integer(labels(colnames(df.AEMB))[-CpG_b])

# Only get the Methylation Data
cat("Grepping only the CpG-data from both dataframes...\n")
df.AEMP.NEW = df.AEMP[,CpG_p]
df.AEMB.NEW = df.AEMB[,CpG_b]

# Filter Probes that are not in the AEM-dataset. Here we get only
# the clinical + QC'ed methylation data of the selected chromosomes.
# Probe_subset_p = Probe_subset[Probe_subset %in% colnames(df.AEMP)]
# Probe_subset_b = Probe_subset[Probe_subset %in% colnames(df.AEMB)]
# df.AEMP.NEW = df.AEMP[,c(colnames(df.AEMP)[clinic_p],Probe_subset_p)]
# df.AEMB.NEW = df.AEMB[,c(colnames(df.AEMB)[clinic_b],Probe_subset_b)]

###-------------------------------------------------------------------------------------------------
### TRANSPOSING NEW DATASETS

cat("Transposing the new dataframes...\n")
df.AEMP.NEW.t = as.data.frame(t(df.AEMP.NEW))
df.AEMB.NEW.t = as.data.frame(t(df.AEMB.NEW))

# Creating CpG Key
# we first create an empty key to which we bind the dataframe, this way it will be at the front,
# instead of at the back
df.AEMP.NEW.t = cbind(IlmnID = NA, df.AEMP.NEW.t)
df.AEMB.NEW.t = cbind(IlmnID = NA, df.AEMB.NEW.t)
df.AEMP.NEW.t$IlmnID <- rownames(df.AEMP.NEW.t)
df.AEMB.NEW.t$IlmnID <- rownames(df.AEMB.NEW.t)

# Select variables from annotation file
myvars <- c("CHR", "MAPINFO", "MAPINFO", "IlmnID")
mtx.HD450K.selection <- mtx.HD450K[myvars]
colnames(mtx.HD450K.selection) <- c("#chr", "start", "end", "IlmnID")

mtx.HD450K.selection$end <- (mtx.HD450K.selection[,3] + 1)

###-------------------------------------------------------------------------------------------------
### MERGING NEW DATASET
cat("Merging the new dataframe with the annotations...\n")
df.AEMP.NEW.t.mtx = merge(mtx.HD450K.selection, df.AEMP.NEW.t, by = "IlmnID", all = FALSE)
df.AEMB.NEW.t.mtx = merge(mtx.HD450K.selection, df.AEMB.NEW.t, by = "IlmnID", all = FALSE)

# Rename column - depends on 'plyr' package
df.AEMP.fastQTL.temp = rename(df.AEMP.NEW.t.mtx, c("IlmnID" = "CpG"))
df.AEMB.fastQTL.temp = rename(df.AEMB.NEW.t.mtx, c("IlmnID" = "CpG"))

# Reorder data
columnsP = seq(5, 492, by = 1)
columnsB = seq(5, 96, by = 1)
df.AEMP.fastQTL.temp2 = df.AEMP.fastQTL.temp[,c(2,1,3,4,columnsP)]
df.AEMB.fastQTL.temp2 = df.AEMB.fastQTL.temp[,c(2,1,3,4,columnsB)]

# Add a leading zero to chromosome 1-9
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "1"] <- "01"
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "2"] <- "02"
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "3"] <- "03"
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "4"] <- "04"
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "5"] <- "05"
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "6"] <- "06"
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "7"] <- "07"
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "8"] <- "08"
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "9"] <- "09"
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "X"] <- "0X"
df.AEMP.fastQTL.temp2$`#chr`[df.AEMP.fastQTL.temp2$`#chr` == "Y"] <- "0Y"

df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "1"] <- "01"
df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "2"] <- "02"
df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "3"] <- "03"
df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "4"] <- "04"
df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "5"] <- "05"
df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "6"] <- "06"
df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "7"] <- "07"
df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "8"] <- "08"
df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "9"] <- "09"
df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "X"] <- "0X"
df.AEMB.fastQTL.temp2$`#chr`[df.AEMB.fastQTL.temp2$`#chr` == "Y"] <- "0Y"

# Sort the data
df.AEMP.fastQTL.temp3 <- as.data.frame(setorder(df.AEMP.fastQTL.temp2, start))
df.AEMP.fastQTL.temp4 <- as.data.frame(setorder(df.AEMP.fastQTL.temp3, `#chr`))
df.AEMB.fastQTL.temp3 <- as.data.frame(setorder(df.AEMB.fastQTL.temp2, start))
df.AEMB.fastQTL.temp4 <- as.data.frame(setorder(df.AEMB.fastQTL.temp3, `#chr`))

# Reorder columns
df.AEMP.fastQTL <- df.AEMP.fastQTL.temp4[ , c(1,3,4,2,5:492)]
df.AEMB.fastQTL <- df.AEMB.fastQTL.temp4[ , c(1,3,4,2,5:96)]


###-------------------------------------------------------------------------------------------------
### SAVE THE DATA
cat("Saving data, while removing intermediate files...\n")

df.AEMP = df.AEMP.original
df.AEMB = df.AEMB.original

rm(df.AEMB.NEW.t.mtx, df.AEMP.NEW.t.mtx, 
   df.AEMB.NEW.t, df.AEMP.NEW.t, 
   df.AEMB.NEW, df.AEMP.NEW,
   df.AEMB.original, df.AEMP.original,
   df.AEMP.fastQTL.temp, df.AEMB.fastQTL.temp,
   df.AEMP.fastQTL.temp2, df.AEMB.fastQTL.temp2,
   df.AEMP.fastQTL.temp3, df.AEMB.fastQTL.temp3,
   df.AEMP.fastQTL.temp4, df.AEMB.fastQTL.temp4)

save.image(paste0(ROOT_loc,"/",Today,"_Create_mQTL_data.RData"))

###-------------------------------------------------------------------------------------------------
### SAVE THE DATA
cat("Writing final data to a text-file...\n")

write.table(df.AEMP.fastQTL, paste0(OUT_loc,"/aems_plaques_450k_QC_443872.bed"), 
            quote = FALSE, sep = "\t", dec = ".", row.names = FALSE)
write.table(df.AEMB.fastQTL, paste0(OUT_loc,"/aems_blood_450k_QC_443872.bed"), 
            quote = FALSE, sep = "\t", dec = ".", row.names = FALSE)

###-------------------------------------------------------------------------------------------------
### DEFINE THE COLORS

###  UtrechtSciencePark Colours Scheme
###
### Website to convert HEX to RGB: http://hex.colorrrs.com.
### For some functions you should divide these numbers by 255.
###
###  Color                         HEX       RGB                 CHR       MAF/INFO
###-------------------------------------------------------------------------------------------------
###	yellow                        #FBB820 (251,184,32)	=>	1 		or 1.0 > INFO
###	gold                          #F59D10 (245,157,16)	=>	2		
###	salmon                        #E55738 (229,87,56) 	=>	3 		or 0.05 < MAF < 0.2 or 0.4 < INFO < 0.6
###	darkpink                      #DB003F ((219,0,63)		=>	4		
###	lightpink                     #E35493 (227,84,147)	=>	5 		or 0.8 < INFO < 1.0
###	pink                          #D5267B (213,38,123)	=>	6		
###	hardpink                      #CC0071 (204,0,113)		=>	7		
###	lightpurple                   #A8448A (168,68,138)	=>	8		
###	purple                        #9A3480 (154,52,128)	=>	9		
###	lavendel                      #8D5B9A (141,91,154)	=>	10		
###	bluepurple                    #705296 (112,82,150)	=>	11		
###	purpleblue                    #686AA9 (104,106,169)	=>	12		
###	lightpurpleblue               #6173AD (97,115,173)	=>	13		
###	seablue                       #4C81BF (76,129,191)	=>	14		
###	skyblue                       #2F8BC9 (47,139,201)	=>	15		
###	azurblue                      #1290D9 (18,144,217)	=>	16		 or 0.01 < MAF < 0.05 or 0.2 < INFO < 0.4
###	lightazurblue                 #1396D8 (19,150,216)	=>	17		
###	greenblue                     #15A6C1 (21,166,193)	=>	18		
###	seaweedgreen                  #5EB17F (94,177,127)	=>	19		
###	yellowgreen                   #86B833 (134,184,51)	=>	20		
###	lightmossgreen                #C5D220 (197,210,32)	=>	21		
###	mossgreen                     #9FC228 (159,194,40)	=>	22		or MAF > 0.20 or 0.6 < INFO < 0.8
###	lightgreen                    #78B113 (120,177,19)	=>	23/X
###	green                         #49A01D (73,160,29)		=>	24/Y
###	grey                          #595A5C (89,90,92)		=>	25/XY	or MAF < 0.01 or 0.0 < INFO < 0.2
###	lightgrey                     #A2A3A4	(162,163,164)	=> 	26/MT
###-------------------------------------------------------------------------------------------------

