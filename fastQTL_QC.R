#!/hpc/local/CentOS7/dhl_ec/software/R-3.2.3/bin/Rscript --vanilla

# Alternative shebang for local Mac OS X: "#!/usr/local/bin/R/bin/Rscript"
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                       fastQTL RESULTS QUALITY CONTROL & PARSER v1_20160303
\n
* Version: fastQTL_QC_v1_20160303
* Last edit: 2016-03-03
* Created by: Sander W. van der Laan | s.w.vanderlaan-2@umcutrecht.nl
\n
* Description:  Results parsing and quality control from fastQTL. The script should be usuable on 
both any Linux distribution with R 3+ installed, Mac OS X and Windows.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

# usage: ./fastQTL_QC.R -p projectdir -r resultfile -o outputdir -t resulttype -q qtltype -a annotfile -j genstatsfile [OPTIONAL: -v verbose (DEFAULT) -q quiet]
#        ./fastQTL_QC.R --projectdir projectdir --resultsfile resultfile --outputdir outputdir --resulttype resulttype --qtltype qtltype --annotfile annotfile --genstats genestatfile [OPTIONAL: --verbose verbose (DEFAULT) -quiet quiet]

# using optparse-library
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
# vignette: http://www.icesi.edu.co/CRAN/web/packages/optparse/vignettes/optparse.pdf
# don't say "Loading required package: optparse"
#suppressPackageStartupMessages(require(optparse))

cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("\nClearing the board, and checking availability of prerequisite packages. We will install packages 
if needed.\n\n")
### CLEAR THE BOARD
rm(list=ls())
### The part of installing (and loading) packages via Rscript doesn't properly work.
### FUNCTION TO INSTALL PACKAGES
install.packages.auto <- function(x) { 
  x <- as.character(substitute(x)) 
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else { 
    # Update installed packages - this may mean a full upgrade of R, which in turn
    # may not be warrented. 
    #update.packages(ask = FALSE) 
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
  }
  if(isTRUE(x %in% .packages(all.available = TRUE))) { 
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
### INSTALL PACKAGES WE NEED
install.packages.auto("optparse")
install.packages.auto("tools")
install.packages.auto("qvalue") # Needed for multiple-testing correction

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###  UtrechtSciencePark Colours Scheme
###
### Website to convert HEX to RGB: http://hex.colorrrs.com.
### For some functions you should divide these numbers by 255.
### 
###	No.	Color				HEX		RGB					CHR		MAF/INFO
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###	1	yellow					    #FBB820 (251,184,32)	=>	1 		or 1.0 > INFO
###	2	gold				        #F59D10 (245,157,16)	=>	2		
###	3	salmon					    #E55738 (229,87,56) 	=>	3 		or 0.05 < MAF < 0.2 or 0.4 < INFO < 0.6
###	4	darkpink				    #DB003F ((219,0,63)		=>	4		
###	5	lightpink				    #E35493 (227,84,147)	=>	5 		or 0.8 < INFO < 1.0
###	6	pink					      #D5267B (213,38,123)	=>	6		
###	7	hardpink				    #CC0071 (204,0,113)		=>	7		
###	8	lightpurple			    #A8448A (168,68,138)	=>	8		
###	9	purple				      #9A3480 (154,52,128)	=>	9		
###	10	lavendel			    #8D5B9A (141,91,154)	=>	10		
###	11	bluepurple			  #705296 (112,82,150)	=>	11		
###	12	purpleblue			  #686AA9 (104,106,169)	=>	12		
###	13	lightpurpleblue		#6173AD (97,115,173/101,120,180)	=>	13		
###	14	seablue				    #4C81BF (76,129,191)	=>	14		
###	15	skyblue				    #2F8BC9 (47,139,201)	=>	15		
###	16	azurblue			    #1290D9 (18,144,217)	=>	16		 or 0.01 < MAF < 0.05 or 0.2 < INFO < 0.4
###	17	lightazurblue		  #1396D8 (19,150,216)	=>	17		
###	18	greenblue			    #15A6C1 (21,166,193)	=>	18		
###	19	seaweedgreen	  	#5EB17F (94,177,127)	=>	19		
###	20	yellowgreen			  #86B833 (134,184,51)	=>	20		
###	21	lightmossgreen		#C5D220 (197,210,32)	=>	21		
###	22	mossgreen			    #9FC228 (159,194,40)	=>	22		or MAF > 0.20 or 0.6 < INFO < 0.8
###	23	lightgreen			  #78B113 (120,177,19)	=>	23/X
###	24	green				      #49A01D (73,160,29)		=>	24/Y
###	25	grey				      #595A5C (89,90,92)		=>	25/XY	or MAF < 0.01 or 0.0 < INFO < 0.2
###	26	lightgrey			    #A2A3A4	(162,163,164)	=> 	26/MT
uithof_color=c("#FBB820","#F59D10","#E55738","#DB003F","#E35493","#D5267B",
               "#CC0071","#A8448A","#9A3480","#8D5B9A","#705296","#686AA9",
               "#6173AD","#4C81BF","#2F8BC9","#1290D9","#1396D8","#15A6C1",
               "#5EB17F","#86B833","#C5D220","#9FC228","#78B113","#49A01D",
               "#595A5C","#A2A3A4")
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
### OPTION LISTING
option_list = list(
  make_option(c("-p", "--projectdir"), action="store", default=NA, type='character',
              help="Path to the project directory."),
  make_option(c("-r", "--resultfile"), action="store", default=NA, type='character',
              help="Path to the results directory, relative to the project directory."),
  make_option(c("-t", "--resulttype"), action="store", default=NA, type='character',
              help="The result type, either [NOM/PERM] for nominal or permutation results, respectively."),
  make_option(c("-q", "--qtltype"), action="store", default=NA, type='character',
              help="The quantitative trait locus (QTL) analysis type , either [EQTL/MQTL] for expression or methylation QTL analysis, respectively."),
  make_option(c("-o", "--outputdir"), action="store", default=NA, type='character',
              help="Path to the output directory."),
  make_option(c("-a", "--annotfile"), action="store", default=NA, type='character',
              help="Path to the annotation file."),
  make_option(c("-j", "--genstats"), action="store", default=NA, type='character',
              help="Path to the summary statistics of the genotypes."),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Should the program print extra stuff out? [default %default]"),
  make_option(c("-s", "--silent"), action="store_false", dest="verbose",
              help="Make the program not be verbose.")
  #make_option(c("-c", "--cvar"), action="store", default="this is c",
  #            help="a variable named c, with a default [default %default]")  
)
opt = parse_args(OptionParser(option_list=option_list))

### FOR LOCAL DEBUGGING
#opt$projectdir="/Users/sanderwvanderlaan/PLINK/analyses/ctmm"
#opt$outputdir="/Users/sanderwvanderlaan/PLINK/analyses/ctmm/" 
#opt$resulttype="PERM"
#opt$qtltype="EQTL"
#opt$resultfile="/Users/sanderwvanderlaan/PLINK/analyses/ctmm/eqtl/rs12539895_CAD/ctmm_QC_eqtlperm_rs12539895_CAD.txt.gz"
#opt$annotfile="/Users/sanderwvanderlaan/PLINK/_CTMM_ORIGINALS/CTMMHumanHT12v4r2_15002873B/annotation_ctmm_all.csv" 
#opt$genstats="/Users/sanderwvanderlaan/PLINK/analyses/ctmm/eqtl/rs12539895_CAD/ctmm_1kGp3GoNL5_QC_rs12539895_CAD.stats"


if (opt$verbose) {
  # if (opt$verbose) {
  # you can use either the long or short name
  # so opt$a and opt$avar are the same.
  # show the user what the variables are
  cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat("Checking the settings.")
  cat("\nThe project directory....................: ")
  cat(opt$projectdir)
  cat("\n\nThe results file.........................: ")
  cat(opt$resultfile)
  cat("\n\nThe output directory.....................: ")
  cat(opt$outputdir)
  cat("\n\nThe annotation file......................: ")
  cat(opt$annotfile)
  cat("\n\nThe results type.........................: ")
  cat(opt$resulttype)
  cat("\n\nThe QTL analysis type....................: ")
  cat(opt$qtltype)
  cat("\n\nThe variant summary statistics...........: ")
  cat(opt$genstats)
  cat("\n\n")
  
}
cat("\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("Wow. We are finally starting \"fastQTL Results Quality Control & Parser v1_20160303\". ")
#--------------------------------------------------------------------------
### START OF THE PROGRAM
# main point of program is here, do this whether or not "verbose" is set
if(!is.na(opt$projectdir) & !is.na(opt$resultfile) & !is.na(opt$outputdir) & !is.na(opt$annotfile) & !is.na(opt$resulttype) & !is.na(opt$qtltype) & !is.na(opt$genstats)) {
  cat(paste("We are going to \nmake some graphs for quality control of you fastQTL analysis. \nThe parsed results from.....: '",file_path_sans_ext(basename(opt$resultfile), compression = TRUE),"'\nand will be outputed in.....: '", opt$outputdir, "'.\n",sep=''))
  cat("\n\n")

  #--------------------------------------------------------------------------
  ### GENERAL SETUP
  Today=format(as.Date(as.POSIXlt(Sys.time())), "%Y%m%d")
  cat(paste("Today's date is: ", Today, ".\n", sep = ''))
  
  #--------------------------------------------------------------------------
  #### DEFINE THE LOCATIONS OF DATA
  ROOT_loc = opt$projectdir # argument 1
  OUT_loc = opt$outputdir # argument 4
  
  #--------------------------------------------------------------------------
  ### LOADING ANNOTATION AND RESULTS FILES DEPENDING ON RESULT TYPE
  cat("Loading annotations...\n")
  ### Location of is set by 'opt$annotfile' # argument 5
  ### The type of the analysis will determine what to load 'opt$qtltype' # argument 4
  CTMMHumanHT12v4r2 = read.csv(opt$annotfile, head = TRUE, stringsAsFactors = FALSE)
  colnames(CTMMHumanHT12v4r2) = c("EntrezID", "ProbeID", "ArrayID", 
                                  "GeneName", "GeneInfo",
                                  "Chr", "GeneTxStart", "GeneTxEnd")

  cat("Loading variant statistics...\n")
  VARIANTSTATS.RAW = read.table(opt$genstats, head = TRUE, stringsAsFactors = FALSE)
  
  # calculate MAC
  VARIANTSTATS.RAW$MAC <- (VARIANTSTATS.RAW[,19]*VARIANTSTATS.RAW[,18]*2)
  
  # calculate caf
  VARIANTSTATS.RAW$CAF <- (((2*VARIANTSTATS.RAW[,16])+VARIANTSTATS.RAW[,15])/(VARIANTSTATS.RAW[,18]*2))
  
  # make imputation column
  VARIANTSTATS.RAW$Imputation <- ifelse(VARIANTSTATS.RAW$alternate_ids == "---", 
                                        c("imputed"), c("genotyped")) 
  
  # Select the columns we need
  VARIANTSTATS = VARIANTSTATS.RAW[,c(2,3,4,5,6, # chr bp
                                     19,        # maf
                                     23,        # mac, column 23
                                     24,        # caf, column 24
                                     8,9,21,18, # imputation quality, HWE and N
                                     25)]       # imputation, column 25
  
  # Change the column names
  colnames(VARIANTSTATS) = c("VARIANT", "Chr", "BP", "OtherAlleleA", "CodedAlleleA", 
                             "MAF", "MAC", "CAF", 
                             "AvgMAxPostCall", "Info", "HWE", "N", "Imputation")
  
  ### Loading nominal results
  if(opt$resulttype == "NOM") { # argument 3
    cat("\nLoading data from 'nominal pass'...\n")
    RESULTS = read.table(opt$resultfile, head = FALSE, stringsAsFactors = FALSE)
    colnames(RESULTS) = c("ProbeID", "VARIANT", "Distance_VARIANT_ProbeID", "Nominal_P", "Beta")
    
    #--------------------------------------------------------------------------
    ### PLOTTING NOMINAL RESULTS
    cat("\nPlotting results...\n") 
    ## To check that the beta approximated permutation p-values are well estimated.
    pdf(paste0(opt$outputdir, # map to the output directory
               Today,"_", # add in Today's date
               file_path_sans_ext(basename(opt$resultfile), compression = TRUE), # get the basename file without the extension and any compression extensions
               "_histogram_nominal_beta.pdf"), onefile = TRUE)
    hist(RESULTS$Beta, 
         breaks = 10000,
         xlab="Effect size", ylab="Distribution", 
         main="Overall distribution of effect size", 
         col = "#1290D9")
    abline(v = mean(RESULTS$Beta), col="#E55738")
    abline(v = (mean(RESULTS$Beta)-4*sd(RESULTS$Beta)), col="#E55738", lty = 2)
    abline(v = (mean(RESULTS$Beta)+4*sd(RESULTS$Beta)), col="#E55738", lty = 2)
    dev.off()
    
  }
  ### Loading permutation results
  else if (opt$resulttype == "PERM") {
    cat("\nLoading data from 'permutation pass'...\n")
    RESULTS = read.table(opt$resultfile, head = FALSE, stringsAsFactors = FALSE)
    colnames(RESULTS) = c("ProbeID", "NVariants", "MLE_Beta_shape1", "MLE_Beta_shape2", "Dummy", 
                          "VARIANT", "Distance_VARIANT_ProbeID", "Nominal_P", "Beta", "Perm_P", "Approx_Perm_P")
    
    #--------------------------------------------------------------------------
    ### PLOTTING NOMINAL RESULTS
    pdf(paste0(opt$outputdir,Today,"_",file_path_sans_ext(basename(opt$resultfile), compression = TRUE),"_comparing_permutation_pvalues.pdf"), onefile = TRUE)
    plot(RESULTS$Perm_P, RESULTS$Approx_Perm_P, 
         xlab="Direct method", ylab="Beta approximation", 
         main="Comparing permuted p-values", bty = "n", 
         pch = 20, col = "#1290D9")
    abline(0, 1, col="#E55738")
    hist(RESULTS$Beta, 
         breaks = 25,
         xlab="Effect size", ylab="Distribution", 
         main="Overall distribution of effect size", 
         #bty = "n", 
         col = "#1290D9"
    )
    abline(v = mean(RESULTS$Beta), col="#E55738")
    abline(v = (mean(RESULTS$Beta)-4*sd(RESULTS$Beta)), col="#E55738", lty = 2)
    abline(v = (mean(RESULTS$Beta)+4*sd(RESULTS$Beta)), col="#E55738", lty = 2)
    dev.off()
    
  }
  else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  }
  
  #--------------------------------------------------------------------------
  ### GET Z-SCORES, SD & SEM
  cat("\nGet Z-scores, sd and sem from p-values...\n")
  ### references:
  ###     - http://stats.stackexchange.com/questions/101136/how-can-i-find-a-z-score-from-a-p-value
  RESULTS$Z = qnorm(RESULTS$Nominal_P)
  
  ### Get standard deviation (SD)
  RESULTS$SD = (RESULTS$Beta-mean(RESULTS$Beta))/RESULTS$Z
  
  ### Get standard error of the mean (SEM)
  RESULTS$SEM = RESULTS$Beta/RESULTS$Z
  
  #--------------------------------------------------------------------------
  #### APPLY MULTIPLE TESTING CORRECTION ###
  cat("\nApplying multiple testing correction methods.\n")
  cat("\n* Conservative correction: Bonferroni correction...\n")
  ### Bonferroni correction - Conservative
  ### references:
  ###     - http://en.wikipedia.org/wiki/Bonferroni_correction
  ###     - https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
  if(opt$resulttype == "NOM")
    RESULTS$Bonferroni = p.adjust(RESULTS$Nominal_P, method="bonferroni")
  if(opt$resulttype == "PERM")
    RESULTS$Bonferroni = p.adjust(RESULTS$Approx_Perm_P, method="bonferroni")
  
  cat("\n* Less conservative correction: Benjamini & Hochberg correction...\n")
  ### Benjamini & Hochberg correction - Less conservative
  ### references:
  ###     - http://en.wikipedia.org/wiki/False_discovery_rate
  ###     - https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
  if(opt$resulttype == "NOM")
    RESULTS$BenjHoch = p.adjust(RESULTS$Nominal_P, method="fdr")
  if(opt$resulttype == "PERM")
    RESULTS$BenjHoch = p.adjust(RESULTS$Approx_Perm_P, method="fdr")
  
  cat("\n* Least conservative correction: Storey & Tibshirani correction...\n")
  ### Storey & Tibshirani correction - Least conservative
  ### references:
  ###     - http://en.wikipedia.org/wiki/False_discovery_rate
  ###     - http://svitsrv25.epfl.ch/R-doc/library/qvalue/html/qvalue.html
  ### Requires a bioconductor package: "qvalue"
  if(opt$resulttype == "NOM")
    RESULTS$Q = qvalue(RESULTS$Nominal_P)$qvalues
  if(opt$resulttype == "PERM")
    RESULTS$Q = qvalue(RESULTS$Approx_Perm_P)$qvalues
  
  #--------------------------------------------------------------------------
  #### ADD IN THE ANNOTATIONS ###
  cat("\nApplying annotations.\n")
  cat("\n* First order based on Benjamini-Hochberg p-values...\n")
  RESULTS.toANNOTATE=RESULTS[order(RESULTS$BenjHoch),]
  
  cat("\n* Now annotate...\n")
  RESULTS.toANNOTATE = cbind(RESULTS.toANNOTATE, CTMMHumanHT12v4r2[match(RESULTS.toANNOTATE[,1], CTMMHumanHT12v4r2$ProbeID ), 
                                                                   c("EntrezID","ArrayID", 
                                                                     "GeneName", "GeneInfo",
                                                                     "Chr", "GeneTxStart", "GeneTxEnd")])
  if(opt$resulttype == "NOM")
    RESULTS.toANNOTATE2 = cbind(RESULTS.toANNOTATE, VARIANTSTATS[match(RESULTS.toANNOTATE[,2], VARIANTSTATS$VARIANT ),])
  if(opt$resulttype == "PERM")
    RESULTS.toANNOTATE2 = cbind(RESULTS.toANNOTATE, VARIANTSTATS[match(RESULTS.toANNOTATE[,6], VARIANTSTATS$VARIANT ),])
  
  ### COLUMNS NOMINAL RESULTS EXAMPLE
  ### NO.	COLUMN-NAME             TYPE  Example
  ### 1.  ProbeID                 : chr  "ILMN_3268799" "ILMN_3268799" "ILMN_3268799" "ILMN_1724504" ...
  ### 2.  VARIANT                 : chr  "rs188294694:99777581:G:A" "14:99807802:C:T" "14:99831389:G:A" "rs1257265" ...  
  ### 3.  Distance_VARIANT_ProbeID: int  -200023 -169802 -146215 1588 7728 13559 16794 21548 23018 25077 ...  
  ### 4.  Nominal_P               : num  1.36e-16 1.36e-16 1.36e-16 3.32e-14 4.74e-14 ...  
  ### 5.  Beta                    : num  -1.3842 -1.3842 -1.3842 0.0966 0.0959 ...  
  ### 6.  Z                       : num  -8.19 -8.19 -8.19 -7.49 -7.45 ...  
  ### 7.  SD                      : num  0.1691 0.1691 0.1691 -0.0129 -0.0129 ...  
  ### 8.  SEM                     : num  0.1691 0.1691 0.1691 -0.0129 -0.0129 ...  
  ### 9.  Bonferroni              : num  1.28e-10 1.28e-10 1.28e-10 3.12e-08 4.46e-08 ...  
  ### 10. BenjHoch                : num  4.26e-11 4.26e-11 4.26e-11 7.36e-10 7.36e-10 ...  
  ### 11. Q                       : num  4.26e-11 4.26e-11 4.26e-11 7.36e-10 7.36e-10 ...  
  ### 12. EntrezID                : int  317762 317762 317762 84193 84193 84193 84193 84193 84193 84193 ...  
  ### 13. ArrayID                 : int  1430092 1430092 1430092 2940706 2940706 2940706 2940706 2940706 2940706 2940706 ...  
  ### 14. GeneName                : chr  "C14orf65" "C14orf65" "C14orf65" "SETD3" ...  
  ### 15. GeneInfo                : chr  "PREDICTED: Homo sapiens misc_RNA (C14orf65), miscRNA." "PREDICTED: Homo sapiens misc_RNA (C14orf65), miscRNA." "PREDICTED: Homo sapiens misc_RNA (C14orf65), miscRNA." "Homo sapiens SET domain containing 3 (SETD3), transcript variant 1, mRNA." ...  
  ### 16. Chr                     : chr  "14" "14" "14" "14" ...  
  ### 17. GeneTxStart             : int  20370725 20821894 17830385 18066400 18535369 18535369 18535369 18535369 19155091 19184405 ...
  ### 18. GeneTxEnd               : int  20455382 20826508 17980131 18067486 19036992 19036992 19036992 19036992 19157295 19185044 ...
  ### 19. VARIANT                 : chr  "rs73684804" "rs185860213" "rs77020184" "rs77451681" ...
  ### 20. Chr                     : int  7 7 7 7 7 7 7 7 7 7 ...
  ### 21. BP                      : int  19370021 19041769 18631905 18944050 19165737 18143868 20715133 18011321 20010778 20367643 ...
  ### 22. OtherAlleleA            : chr  "G" "T" "T" "T" ...
  ### 23. CodedAlleleA            : chr  "T" "G" "C" "C" ...
  ### 24. MAF                     : num  0.01787 0.00522 0.11404 0.06209 0.04491 ...
  ### 25. MAC                     : num  22.3 6.51 142.32 77.49 56.04 ...
  ### 26. CAF                     : num  0.01787 0.00522 0.11404 0.06209 0.04491 ...
  ### 27. AvgMAxPostCall          : num  0.997 0.999 0.985 0.987 0.993 ...
  ### 28. Info                    : num  0.935 0.941 0.944 0.908 0.932 ...
  ### 29. HWE                     : num  1 1 1 0.718 0.338 ...
  ### 30. N                       : int  624 624 624 624 624 624 624 624 624 624 ...
  ### 31. Imputation              : chr  "imputed" "imputed" "imputed" "imputed" ...
  
  ### COLUMNS PERMUTATION RESULTS EXAMPLE
  ### 1.  ProbeID                 : chr  "ILMN_1749930" "ILMN_1739594" "ILMN_1696827" "ILMN_1726420" ...
  ### 2.  NVariants               : int  16252 16252 16252 16252 16252 16252 16252 16252 16252 16252 ...
  ### 3.  MLE_Beta_shape1         : num  0.894 1.022 1.042 1.008 1.045 ...
  ### 4.  MLE_Beta_shape2         : num  359 1391 1136 782 1686 ...
  ### 5.  Dummy                   : num  189 259 249 220 276 ...
  ### 6.  VARIANT                 : chr  "chr1:55815805:D" "rs33973564" "rs77205069" "rs192525260" ...
  ### 7.  Distance_VARIANT_ProbeID: int  1584670 -342113 15598 -931086 1210570 216139 3336388 132711 1318328 272206 ...
  ### 8.  Nominal_P               : num  1.03e-06 9.18e-06 8.85e-06 9.22e-07 4.01e-05 ...
  ### 9.  Beta                    : num  -0.2993 0.0977 -0.0633 -0.4006 -0.3329 ...
  ### 10. Perm_P                  : num  0.0856 0.0417 0.0445 0.0201 0.0999 ...
  ### 11. Approx_Perm_P           : num  0.0507 0.0414 0.0441 0.0172 0.1044 ...
  ### 12. Z                       : num  -4.75 -4.28 -4.29 -4.77 -3.94 ...
  ### 13. SD                      : num  0.032 -0.0572 -0.0196 0.0531 0.047 ...
  ### 14. SEM                     : num  0.0631 -0.0228 0.0148 0.084 0.0844 ...
  ### 15. Bonferroni              : num  1 1 1 0.722 1 ...
  ### 16. BenjHoch                : num  0.532 0.532 0.532 0.532 0.583 ...
  ### 17. Q                       : num  0.238 0.238 0.238 0.238 0.26 ...
  ### 18. EntrezID                : int  55706 26027 25973 7809 63948 54432 1733 115353 127428 200008 ...
  ### 19. ArrayID                 : int  4200601 1050253 3400170 3440435 1070341 1400133 1470164 5570524 7510553 6180367 ...
  ### 20. GeneName                : chr  "TMEM48" "ACOT11" "PARS2" "BSND" ...
  ### 21. GeneInfo                : chr  "Homo sapiens transmembrane protein 48 (TMEM48), mRNA." "Homo sapiens acyl-CoA thioesterase 11 (ACOT11), transcript variant 2, mRNA." "Homo sapiens prolyl-tRNA synthetase 2, mitochondrial (putative) (PARS2), nuclear gene encoding mitochondrial protein, mRNA." "Homo sapiens Bartter syndrome, infantile, with sensorineural deafness (Barttin) (BSND), mRNA." ...
  ### 22. Chr                     : chr  "1" "1" "1" "1" ...
  ### 23. GeneTxStart             : int  20370725 20821894 17830385 18066400 18535369 18535369 18535369 18535369 19155091 19184405 ...
  ### 24. GeneTxEnd               : int  20455382 20826508 17980131 18067486 19036992 19036992 19036992 19036992 19157295 19185044 ...
  ### 25. VARIANT                 : chr  "rs73684804" "rs185860213" "rs77020184" "rs77451681" ...
  ### 26. Chr                     : int  7 7 7 7 7 7 7 7 7 7 ...
  ### 27. BP                      : int  19370021 19041769 18631905 18944050 19165737 18143868 20715133 18011321 20010778 20367643 ...
  ### 28. OtherAlleleA            : chr  "G" "T" "T" "T" ...
  ### 29. CodedAlleleA            : chr  "T" "G" "C" "C" ...
  ### 30. MAF                     : num  0.01787 0.00522 0.11404 0.06209 0.04491 ...
  ### 31. MAC                     : num  22.3 6.51 142.32 77.49 56.04 ...
  ### 32. CAF                     : num  0.01787 0.00522 0.11404 0.06209 0.04491 ...
  ### 33. AvgMAxPostCall          : num  0.997 0.999 0.985 0.987 0.993 ...
  ### 34. Info                    : num  0.935 0.941 0.944 0.908 0.932 ...
  ### 35. HWE                     : num  1 1 1 0.718 0.338 ...
  ### 36. N                       : int  624 624 624 624 624 624 624 624 624 624 ...
  ### 37. Imputation              : chr  "imputed" "imputed" "imputed" "imputed" ...
  
  if(opt$resulttype == "NOM")
    RESULTS.ANNOTATE = RESULTS.toANNOTATE2[,c(1,2,20,21,22,23,24,25,26,29,28,31,30, # Variant information
                                              14,12,3,16,17,18, # Gene information
                                              5,8,4,9,10,11)] # statistics
  if(opt$resulttype == "PERM")
    RESULTS.ANNOTATE = RESULTS.toANNOTATE2[,c(1,6,26,27,28,29,30,31,32,35,34,37,36, # Variant information
                                              20,18,7,22,23,24, # Gene information
                                              9,14,8,10,11,15,16,17)] # statistics
  
  cat("\n* Remove duplicate gene names...\n")
  RESULTS.ANNOTATE[, "GeneName"] = as.character(lapply(RESULTS.toANNOTATE2[,"GeneName"], 
                                                       FUN = function(x){paste(unique(unlist(strsplit(x, split = ";"))), sep="", collapse=";")}))
  
  cat("\n* Correct Colnames...\n")
  if(opt$resulttype == "NOM")
    colnames(RESULTS.ANNOTATE) = c("ProbeID", "VARIANT", "Chr", "BP", "OtherAlleleA", "CodedAlleleA", "MAF", "MAC", "CAF", "HWE", "Info", "Imputation", "N", 
                                   "GeneName", "EntrezID", "Distance_VARIANT_GENE", "Chr", "GeneTxStart", "GeneTxEnd",
                                   "Beta", "SE", "Nominal_P", "Bonferroni","BenjHoch","Q")
  if(opt$resulttype == "PERM")
    colnames(RESULTS.ANNOTATE) = c("ProbeID", "VARIANT", "Chr", "BP", "OtherAlleleA", "CodedAlleleA", "MAF", "MAC", "CAF", "HWE", "Info", "Imputation", "N", 
                                   "GeneName","EntrezID", "Distance_VARIANT_GENE", "Chr", "GeneTxStart", "GeneTxEnd",
                                   "Beta", "SE", "Nominal_P","Perm_P","ApproxPerm_P", "Bonferroni","BenjHoch","Q")
  
  cat("\n* Remove temporary files...\n")
  rm(RESULTS.toANNOTATE, RESULTS.toANNOTATE2)
  
  #--------------------------------------------------------------------------
  ### SAVE NEW DATA ###
  if(opt$resulttype == "NOM")
    #write.table(RESULTS.ANNOTATE[which(RESULTS.ANNOTATE$Q <= 0.05), ], # with filtering on Q-value 
    write.table(RESULTS.ANNOTATE[which(RESULTS.ANNOTATE$Q != "NA"), ], # without filtering on Q-value
                paste0(opt$outputdir, Today,"_", file_path_sans_ext(basename(opt$resultfile), compression = TRUE), 
                       "_nominal.P0_05.txt"), 
                quote = FALSE , row.names = FALSE, col.names = TRUE, sep = "\t")
  if(opt$resulttype == "PERM")
    write.table(RESULTS.ANNOTATE[which(RESULTS.ANNOTATE$Q <= 0.05), ], 
                paste0(opt$outputdir, Today,"_", file_path_sans_ext(basename(opt$resultfile), compression = TRUE), 
                       "_perm.P0_05.txt"), 
                quote = FALSE , row.names = FALSE, col.names = TRUE, sep = "\t")
  
} else {
  cat("You didn't specify all variables:\n
      - --p/projectdir : path to project directory\n
      - --r/resultdir  : path to results directory\n
      - --o/outputdir  : path to output directory\n
      - --t/resulttype : the results type (NOM for nominal; PERM for permutation)\n
      - --q/qtltype    : the QTL analysis type (EQTL for expression QTL; MQTL for methylation QTL)\n
      - --a/annotfile  : path to annotation file of genes\n
      - --j/genstats   : path to summary statistics of variants\n\n", 
      file=stderr()) # print error messages to stderr
}

#--------------------------------------------------------------------------
### CLOSING MESSAGE
cat(paste("All done parsing fastQTL data on",file_path_sans_ext(basename(opt$resultfile), compression = TRUE),".\n"))
cat(paste("\nToday's: ",Today, "\n"))
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")


#--------------------------------------------------------------------------
### SAVE ENVIRONMENT | FOR DEBUGGING
if(opt$resulttype == "NOM")
  save.image(paste0(opt$outputdir,Today,"_",file_path_sans_ext(basename(opt$resultfile), compression = TRUE),"_NOM_DEBUG_FastQTL_analysis.RData"))
if(opt$resulttype == "PERM")
  save.image(paste0(opt$outputdir,Today,"_",file_path_sans_ext(basename(opt$resultfile), compression = TRUE),"_PERM_DEBUG_FastQTL_analysis.RData"))

