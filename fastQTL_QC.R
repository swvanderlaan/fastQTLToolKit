#!/hpc/local/CentOS7/dhl_ec/software/R-3.4.0/bin/Rscript --vanilla

# Alternative shebang for local Mac OS X: "#!/usr/local/bin/Rscript --vanilla"
# Linux version for HPC: #!/hpc/local/CentOS7/dhl_ec/software/R-3.4.0/bin/Rscript --vanilla
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fastQTL RESULTS QUALITY CONTROL & PARSER v2
    \n
    * Version: v2.2.3
    * Last edit: 2018-03-12
    * Created by: Sander W. van der Laan | s.w.vanderlaan-2@umcutrecht.nl
    \n
    * Description:  Results parsing and quality control from fastQTL results using CTMM (eQTL) or 
    Athero-Express (mQTL) data. The script should be usuable on both any Linux distribution with 
    R 3+ installed, Mac OS X and Windows.
    
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

# usage: ./fastQTL_QC.R -p projectdir -r resultfile -o outputdir -t resulttype -q qtltype -a annotfile -j genstatsfile [OPTIONAL: -v verbose (DEFAULT) -q quiet]
#        ./fastQTL_QC.R --projectdir projectdir --resultsfile resultfile --outputdir outputdir --resulttype resulttype --qtltype qtltype --annotfile annotfile --genstats genestatfile [OPTIONAL: --verbose verbose (DEFAULT) -quiet quiet]

#--------------------------------------------------------------------------
cat("\n* Clearing the environment...\n\n")
### CLEAR THE BOARD
rm(list=ls())

cat("\n* Loading function to install packages...\n\n")
### Prerequisite: 'optparse'-library
### * Manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
### * Vignette: http://www.icesi.edu.co/CRAN/web/packages/optparse/vignettes/optparse.pdf

### Don't say "Loading required package: optparse"...
###suppressPackageStartupMessages(require(optparse))
###require(optparse)

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
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE, repos = \"http://cran-mirror.cs.uu.nl/\")", x)))
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

cat("\n* Checking availability of required packages and installing if needed...\n\n")
### INSTALL PACKAGES WE NEED
install.packages.auto("optparse")
install.packages.auto("tools")
install.packages.auto("qvalue") # Needed for multiple-testing correction

cat("\nDone! Required packages installed and loaded.\n\n")

cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

cat("\n* Setting colours...\n\n")
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

### OPTIONLIST | FOR LOCAL DEBUGGING
# opt$projectdir="/Volumes/MyBookStudioII/Backup/PLINK/analyses/test_qtl/mqtl_aems450k1/"
# opt$outputdir="/Volumes/MyBookStudioII/Backup/PLINK/analyses/test_qtl/mqtl_aems450k1/qtl/rs10953541_CAD"
# opt$resulttype="NOM"
# #opt$resulttype="PERM"
# opt$qtltype="MQTL"
# opt$resultfile="/Volumes/MyBookStudioII/Backup/PLINK/analyses/test_qtl/mqtl_aems450k1/qtl/rs10953541_CAD/aegs_QC_qtlnom_rs10953541_excl_DEFAULT.txt.gz"
# #opt$resultfile="/Volumes/MyBookStudioII/Backup/PLINK/analyses/test_qtl/mqtl_aems450k1/qtl/rs10953541_CAD/aegs_QC_qtlperm_rs10953541_excl_DEFAULT.txt.gz"
# opt$annotfile="/Volumes/MyBookStudioII/Backup/PLINK/_AE_Originals/IlluminaMethylation450K.annotation.txt.gz"
# opt$genstats="/Volumes/MyBookStudioII/Backup/PLINK/analyses/test_qtl/mqtl_aems450k1/qtl/rs10953541_CAD/aegs_1kGp3GoNL5_QC_rs10953541_excl_DEFAULT.stats"
# 
# opt$projectdir="/Volumes/MyBookStudioII/Backup/PLINK/analyses/test_qtl/eqtl_ctmm/"
# opt$outputdir="/Volumes/MyBookStudioII/Backup/PLINK/analyses/test_qtl/eqtl_ctmm/qtl/rs10953541_CAD"
# opt$resulttype="NOM"
# #opt$resulttype="PERM"
# opt$qtltype="EQTL"
# opt$resultfile="/Volumes/MyBookStudioII/Backup/PLINK/analyses/test_qtl/eqtl_ctmm/qtl/rs10953541_CAD/ctmm_QC_qtlnom_rs10953541_excl_DEFAULT.txt.gz"
# #opt$resultfile="/Volumes/MyBookStudioII/Backup/PLINK/analyses/test_qtl/eqtl_ctmm/qtl/rs10953541_CAD/ctmm_QC_qtlperm_rs10953541_excl_DEFAULT.txt.gz"
# opt$annotfile="/Volumes/MyBookStudioII/Backup/PLINK/_CTMM_Originals/CTMMHumanHT12v4r2_15002873B/annotation_ctmm_all.txt"
# opt$genstats="/Volumes/MyBookStudioII/Backup/PLINK/analyses/test_qtl/eqtl_ctmm/qtl/rs10953541_CAD/ctmm_1kGp3GoNL5_QC_rs10953541_excl_DEFAULT.stats"
# 
### OPTIONLIST | FOR LOCAL DEBUGGING

if (opt$verbose) {
  # You can use either the long or short name; so opt$a and opt$avar are the same.
  # Show the user what the variables are.
  cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat("* Checking the settings as given through the flags.")
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
  cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  cat("\n\n")
}
cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
cat("Wow, we are all set. Starting \"fastQTL Results Quality Control & Parser\".")
#--------------------------------------------------------------------------
### START OF THE PROGRAM
# main point of program is here, do this whether or not "verbose" is set
if(!is.na(opt$projectdir) & !is.na(opt$resultfile) & !is.na(opt$outputdir) & !is.na(opt$annotfile) & !is.na(opt$resulttype) & !is.na(opt$qtltype) & !is.na(opt$genstats)) {
  cat(paste("\nWe are going to make some graphs for quality control of you fastQTL analysis. \n\nAnalysing these results...............: '",file_path_sans_ext(basename(opt$resultfile), compression = TRUE),"'\nParsed results will be saved here.....: '", opt$outputdir, "'.\n",sep=''))
  
  #--------------------------------------------------------------------------
  ### GENERAL SETUP
  Today=format(as.Date(as.POSIXlt(Sys.time())), "%Y%m%d")
  cat(paste("\nToday's date is: ", Today, ".\n", sep = ''))
  
  #--------------------------------------------------------------------------
  #### DEFINE THE LOCATIONS OF DATA
  ROOT_loc = opt$projectdir # argument 1
  OUT_loc = opt$outputdir # argument 4
  
  cat("\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
  #--------------------------------------------------------------------------
  ### LOADING ANNOTATION AND RESULTS FILES DEPENDING ON RESULT TYPE
  cat("\nLoading annotations...\n")
  ### Location of is set by 'opt$annotfile' # argument 5
  ### The type of the analysis will determine what to load 'opt$qtltype' # argument 4
  if(opt$qtltype == "EQTL") { 
    cat ("\n...for a CTMM based eQTL analysis in monocytes...\n")
    ANNOTATIONSFILE = read.csv(opt$annotfile, head = TRUE, stringsAsFactors = FALSE, sep = ",")
    colnames(ANNOTATIONSFILE) = c("EntrezID", "ProbeID", "ArrayID", 
                                  "GeneName", "GeneInfo","Chr", "GeneTxStart", "GeneTxEnd")
  } else if (opt$qtltype == "MQTL") {
    cat ("\n...for an Athero-Express based MQTL analysis...\n")
    ANNOTATIONSFILE = read.table(opt$annotfile, head = TRUE, stringsAsFactors = FALSE, sep = ",", na.strings = "")
    
    colnames(ANNOTATIONSFILE) = c("IlmnID", "ProbeID", 
                                  "AddressA_ID", "AlleleA_ProbeSeq", "AddressB_ID", "AlleleB_ProbeSeq", 
                                  "Infinium_Design_Type", "Next_Base", "Color_Channel", "Forward_Sequence", 
                                  "Genome_Build", "CHR", "MAPINFO", "SourceSeq", "Chromosome_36", "Coordinate_36", "Strand", 
                                  "Probe_SNPs", "Probe_SNPs_10", "Random_Loci", "Methyl27_Loci", 
                                  "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island", 
                                  "Phantom", "DMR", "Enhancer", "HMM_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DHS", 
                                  "UCSC_RefGene_Dist")
  } else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  }
  
  cat("\nLoading variant statistics...\n")
  VARIANTSTATS.RAW = read.table(opt$genstats, head = TRUE, stringsAsFactors = FALSE)
  cat("\n* calculating 'minor allele count' (MAC)...")
  # calculate MAC
  VARIANTSTATS.RAW$MAC <- (VARIANTSTATS.RAW[,19]*VARIANTSTATS.RAW[,18]*2)
  
  cat("\n* calculating 'coded allele frequency' (CAF)...")
  # calculate caf
  VARIANTSTATS.RAW$CAF <- (((2*VARIANTSTATS.RAW[,16])+VARIANTSTATS.RAW[,15])/(VARIANTSTATS.RAW[,18]*2))
  
  cat("\n* determining which variants are solely 'imputed'...")
  # make imputation column
  VARIANTSTATS.RAW$Imputation <- ifelse(VARIANTSTATS.RAW$alternate_ids == "---", 
                                        c("imputed"), c("genotyped")) 
  
  cat("\n* selecting required variant statistics data...")
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
  
  ### Loading *nominal* results
  if(opt$resulttype == "NOM") { # argument 3
    cat("\n\nLoading data from 'nominal pass'...\n")
    RESULTS = read.table(opt$resultfile, head = FALSE, stringsAsFactors = FALSE)
    colnames(RESULTS) = c("ProbeID", "VARIANT", "Distance_VARIANT_ProbeID", "Nominal_P", "Beta")
    
    #--------------------------------------------------------------------------
    ### PLOTTING NOMINAL RESULTS
    cat("\nPlotting results...\n") 
    ## To check that the beta approximated permutation p-values are well estimated.
    pdf(paste0(opt$outputdir, "/",# map to the output directory
               ###Today,"_", # add in Today's date -- removed as it causes issues in downstream projects when its the 'next day'
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
    
  } else if (opt$resulttype == "PERM") { ### Loading *permutation* results 
    cat("\nLoading data from 'permutation pass'...\n")
    RESULTS = read.table(opt$resultfile, head = FALSE, stringsAsFactors = FALSE)
    colnames(RESULTS) = c("ProbeID", "NVariants", "MLE_Beta_shape1", "MLE_Beta_shape2", "Dummy", 
                          "VARIANT", "Distance_VARIANT_ProbeID", "Nominal_P", "Beta", "Perm_P", "Approx_Perm_P")
    
    #--------------------------------------------------------------------------
    ### PLOTTING PERMUTATION RESULTS
    pdf(paste0(opt$outputdir, "/",# map to the output directory
               ###Today,"_", # add in Today's date -- removed as it causes issues in downstream projects when its the 'next day'
               file_path_sans_ext(basename(opt$resultfile), compression = TRUE), # get the basename file without the extension and any compression extensions
               "_comparing_permutation_pvalues.pdf"), onefile = TRUE)
    
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
    
  } else {
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
  if(opt$resulttype == "NOM") {
    RESULTS$Bonferroni = p.adjust(RESULTS$Nominal_P, method="bonferroni")
  } else if(opt$resulttype == "PERM"){
    RESULTS$Bonferroni = p.adjust(RESULTS$Approx_Perm_P, method="bonferroni")
  } else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  }
  
  cat("\n* Less conservative correction: Benjamini & Hochberg correction...\n")
  ### Benjamini & Hochberg correction - Less conservative
  ### references:
  ###     - http://en.wikipedia.org/wiki/False_discovery_rate
  ###     - https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
  if(opt$resulttype == "NOM") {
    RESULTS$BenjHoch = p.adjust(RESULTS$Nominal_P, method="fdr")
  } else if(opt$resulttype == "PERM") {
    RESULTS$BenjHoch = p.adjust(RESULTS$Approx_Perm_P, method="fdr")
  } else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  }
  
  cat("\n* Least conservative correction: Storey & Tibshirani correction...\n")
  ### Storey & Tibshirani correction - Least conservative
  ### references:
  ###     - http://en.wikipedia.org/wiki/False_discovery_rate
  ###     - http://svitsrv25.epfl.ch/R-doc/library/qvalue/html/qvalue.html
  ### Requires a bioconductor package: "qvalue"
  if(opt$resulttype == "NOM") {
    RESULTS$Q = qvalue(RESULTS$Nominal_P)$qvalues
  } else if(opt$resulttype == "PERM") {
    RESULTS$Q = qvalue(RESULTS$Approx_Perm_P)$qvalues
  } else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  }
  
  #--------------------------------------------------------------------------
  #### ADD IN THE ANNOTATIONS ###
  cat("\nApplying annotations.\n")
  cat("\n* First order based on Benjamini-Hochberg p-values...\n")
  RESULTS.toANNOTATE=RESULTS[order(RESULTS$BenjHoch),]
  
  cat("\n* Now annotating...\n")
  if(opt$qtltype == "EQTL") { 
    cat ("\n...the results of a CTMM based eQTL analysis in monocytes.\n")
    RESULTS.toANNOTATE = cbind(RESULTS.toANNOTATE, ANNOTATIONSFILE[match(RESULTS.toANNOTATE[,1], ANNOTATIONSFILE$ProbeID ), 
                                                                   c("EntrezID","ArrayID", 
                                                                     "GeneName", "GeneInfo",
                                                                     "Chr", "GeneTxStart", "GeneTxEnd")])
    
  } else if (opt$qtltype == "MQTL") {
    cat ("\n...the results of an Athero-Express based MQTL analysis.\n")
    RESULTS.toANNOTATE = cbind(RESULTS.toANNOTATE, ANNOTATIONSFILE[match(RESULTS.toANNOTATE[,1], ANNOTATIONSFILE$ProbeID ), 
                                                                   c("IlmnID", "ProbeID", 
                                                                     "AddressA_ID", "AlleleA_ProbeSeq", "AddressB_ID", "AlleleB_ProbeSeq", 
                                                                     "Infinium_Design_Type", "Next_Base", "Color_Channel", "Forward_Sequence", 
                                                                     "Genome_Build", "CHR", "MAPINFO", "SourceSeq", "Chromosome_36", "Coordinate_36", "Strand", 
                                                                     "Probe_SNPs", "Probe_SNPs_10", "Random_Loci", "Methyl27_Loci", 
                                                                     "UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island", 
                                                                     "Phantom", "DMR", "Enhancer", "HMM_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DHS", 
                                                                     "UCSC_RefGene_Dist")])
  } else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  }
  
  cat("\n* Merging results, genetic stats, and annotations...\n")
  if(opt$resulttype == "NOM") {
    RESULTS.toANNOTATE2 = cbind(RESULTS.toANNOTATE, VARIANTSTATS[match(RESULTS.toANNOTATE[,2], VARIANTSTATS$VARIANT ),])
  } else if(opt$resulttype == "PERM") {
    RESULTS.toANNOTATE2 = cbind(RESULTS.toANNOTATE, VARIANTSTATS[match(RESULTS.toANNOTATE[,6], VARIANTSTATS$VARIANT ),])
  } else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  }
  
  if(opt$qtltype == "EQTL") { 
    cat ("\n* Parsing annotated results for a CTMM eQTL analysis in monocytes...\n")
    if(opt$resulttype == "NOM") {
      cat ("\n--- nominal results ---\n")
      RESULTS.ANNOTATE = RESULTS.toANNOTATE2[,c(1,2,20,21,22,23,24,25,26,29,28,31,30, # Variant information
                                                14,12,3,16,17,18, # Gene information
                                                5,8,4,9,10,11)] # association statistics
    } else if(opt$resulttype == "PERM") {
      cat ("\n--- permuted results ---\n")
      RESULTS.ANNOTATE = RESULTS.toANNOTATE2[,c(1,6,26,27,28,29,30,31,32,35,34,37,36, # Variant information
                                                20,18,7,22,23,24, # Gene information
                                                9,14,8,10,11,15,16,17)] # association statistics
    } else {
      cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
           file=stderr()) # print error messages to stder
    }
    
  } else if(opt$qtltype == "MQTL") {
    cat ("\n* Parsing annotated results for an Athero-Express mQTL analysis...\n")
    if(opt$resulttype == "NOM") {
      cat ("\n--- nominal results ---\n")
      RESULTS.ANNOTATE = RESULTS.toANNOTATE2[,c(1,2,47,48,49,50,51,52,53,56,55,58,57, # Variant information
                                                3,23,24,18, # CpG information
                                                33,34,35,37,38,39,40,41,42,43,44, # CpG associated information
                                                5,8,4,9,10,11)] # association statistics
    } else if(opt$resulttype == "PERM") {
      cat ("\n--- permuted results ---\n")
      RESULTS.ANNOTATE = RESULTS.toANNOTATE2[,c(1,6,53,54,55,56,57,58,59,62,61,64,63, # Variant information
                                                7,29,30,24, # CpG information
                                                39,40,41,43,44,45,46,47,48,49,50, # CpG associated information
                                                9,14,8,10,11,15,16,17)] # association statistics
    } else {
      cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
           file=stderr()) # print error messages to stder
    }
    
  } else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  } 
  
  cat("\n* Remove duplicate gene names...\n")
  if(opt$qtltype == "EQTL") { 
    cat ("\n...for results of a CTMM eQTL analysis in monocytes...\n")
    RESULTS.ANNOTATE[, "GeneName"] = as.character(lapply(RESULTS.toANNOTATE2[,"GeneName"], 
                                                         FUN = function(x){paste(unique(unlist(strsplit(x, split = ";"))), sep="", collapse=";")}))
  } else if(opt$qtltype == "MQTL") {
    cat ("\n...for results of an Athero-Express mQTL analysis...\n")
    RESULTS.ANNOTATE[, "UCSC_RefGene_Name"] = as.character(lapply(RESULTS.toANNOTATE2[,"UCSC_RefGene_Name"], 
                                                                  FUN = function(x){paste(unique(unlist(strsplit(x, split = ";"))), sep="", collapse=";")}))
    RESULTS.ANNOTATE[, "UCSC_RefGene_Accession"] = as.character(lapply(RESULTS.toANNOTATE2[,"UCSC_RefGene_Accession"],
                                                                       FUN = function(x){paste(unique(unlist(strsplit(x, split = ";"))), sep="", collapse=";")}))
    RESULTS.ANNOTATE[, "UCSC_RefGene_Group"] = as.character(lapply(RESULTS.toANNOTATE2[,"UCSC_RefGene_Group"],
                                                                   FUN = function(x){paste(unique(unlist(strsplit(x, split = ";"))), sep="", collapse=";")}))
  } else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  }
  
  cat("\n* Correct Colnames...\n")
  if(opt$qtltype == "EQTL") { 
    cat ("\n...for results of a CTMM eQTL analysis in monocytes...\n")
    if(opt$resulttype == "NOM") {
      colnames(RESULTS.ANNOTATE) = c("ProbeID", "VARIANT", "Chr", "BP", "OtherAlleleA", "CodedAlleleA", "MAF", "MAC", "CAF", "HWE", "Info", "Imputation", "N", 
                                     "GeneName", "EntrezID", "Distance_VARIANT_GENE", "Chr_Gene", "GeneTxStart", "GeneTxEnd",
                                     "Beta", "SE", "Nominal_P", "Bonferroni","BenjHoch","Q")
    } else if(opt$resulttype == "PERM") {
      colnames(RESULTS.ANNOTATE) = c("ProbeID", "VARIANT", "Chr", "BP", "OtherAlleleA", "CodedAlleleA", "MAF", "MAC", "CAF", "HWE", "Info", "Imputation", "N", 
                                     "GeneName","EntrezID", "Distance_VARIANT_GENE", "Chr_Gene", "GeneTxStart", "GeneTxEnd",
                                     "Beta", "SE", "Nominal_P","Perm_P","ApproxPerm_P", "Bonferroni","BenjHoch","Q")
    } else {
      cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
           file=stderr()) # print error messages to stder
    }
    
  } else if(opt$qtltype == "MQTL") {
    cat ("\n...for results of an Athero-Express mQTL analysis...\n")
    if(opt$resulttype == "NOM") {
      colnames(RESULTS.ANNOTATE) = c("ProbeID", "VARIANT", "Chr", "BP", "OtherAlleleA", "CodedAlleleA", "MAF", "MAC", "CAF", "HWE", "Info", "Imputation", "N", 
                                     "Distance_VARIANT_CpG", "Chr_CpG", "BP_CpG",
                                     "ProbeType", "GeneName_UCSC", "AccessionID_UCSC", "GeneGroup_UCSC", 
                                     "CpG_Island_Relation_UCSC", "Phantom", "DMR", "Enhancer", "HMM_Island",
                                     "RegulatoryFeatureName", "RegulatoryFeatureGroup", "DHS",
                                     "Beta", "SE", "Nominal_P", "Bonferroni","BenjHoch","Q")
    } else if(opt$resulttype == "PERM") {
      colnames(RESULTS.ANNOTATE) = c("ProbeID", "VARIANT", "Chr", "BP", "OtherAlleleA", "CodedAlleleA", "MAF", "MAC", "CAF", "HWE", "Info", "Imputation", "N", 
                                     "Distance_VARIANT_CpG", "Chr_CpG", "BP_CpG",
                                     "ProbeType", "GeneName_UCSC", "AccessionID_UCSC", "GeneGroup_UCSC", 
                                     "CpG_Island_Relation_UCSC", "Phantom", "DMR", "Enhancer", "HMM_Island",
                                     "RegulatoryFeatureName", "RegulatoryFeatureGroup", "DHS",
                                     "Beta", "SE", "Nominal_P","Perm_P","ApproxPerm_P", "Bonferroni","BenjHoch","Q")
    } else {
      cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
           file=stderr()) # print error messages to stder
    }
    
  } else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  }
  
  cat("\n* Remove temporary files...\n")
  rm(RESULTS.toANNOTATE, RESULTS.toANNOTATE2)
  
  #--------------------------------------------------------------------------
  ### SAVE NEW DATA ###
  cat("\n* Saving parsed data...\n")
  if(opt$resulttype == "NOM") {
    #write.table(RESULTS.ANNOTATE[which(RESULTS.ANNOTATE$Q <= 0.05), ], # with filtering on Q-value 
    write.table(RESULTS.ANNOTATE[which(RESULTS.ANNOTATE$Q != "NA"), ], # without filtering on Q-value
                #paste0(opt$outputdir, "/", 
                #       ###Today,"_", # add in Today's date -- removed as it causes issues in downstream projects when its the 'next day'
                #       file_path_sans_ext(basename(opt$resultfile), compression = TRUE), 
                #       "_nominal.P0_05.txt"),
                paste0(opt$outputdir, "/", 
                       ###Today,"_", # add in Today's date -- removed as it causes issues in downstream projects when its the 'next day'
                       file_path_sans_ext(basename(opt$resultfile), compression = TRUE), 
                       "_nominal.all.txt"),
                quote = FALSE , row.names = FALSE, col.names = TRUE, sep = ",", na = "NA", dec = ".")
  } else if(opt$resulttype == "PERM") {
    write.table(RESULTS.ANNOTATE[which(RESULTS.ANNOTATE$Q <= 0.05), ], # with filtering on Q-value 
                #write.table(RESULTS.ANNOTATE[which(RESULTS.ANNOTATE$Q != "NA"), ], # without filtering on Q-value
                paste0(opt$outputdir, "/", 
                       ###Today,"_", # add in Today's date -- removed as it causes issues in downstream projects when its the 'next day'
                       file_path_sans_ext(basename(opt$resultfile), compression = TRUE), 
                       "_perm.P0_05.txt"),
                #paste0(opt$outputdir, "/", 
                #        ###Today,"_", # add in Today's date -- removed as it causes issues in downstream projects when its the 'next day' 
                #        file_path_sans_ext(basename(opt$resultfile), compression = TRUE), 
                #       "_perm.all.txt"),
                quote = FALSE , row.names = FALSE, col.names = TRUE, sep = ",", na = "NA", dec = ".")
  } else {
    cat ("\n\n*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please.\n\n", 
         file=stderr()) # print error messages to stder
  }
  
} else {
  cat("*** ERROR *** You didn't specify all variables:\n
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
cat("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")

#--------------------------------------------------------------------------
### SAVE ENVIRONMENT | FOR DEBUGGING
# if(opt$resulttype == "NOM")
#   save.image(paste0(opt$outputdir, "/",Today,"_",file_path_sans_ext(basename(opt$resultfile), compression = TRUE),"_NOM_DEBUG_FastQTL_analysis.RData"))
# if(opt$resulttype == "PERM")
#   save.image(paste0(opt$outputdir, "/",Today,"_",file_path_sans_ext(basename(opt$resultfile), compression = TRUE),"_PERM_DEBUG_FastQTL_analysis.RData"))



###	UtrechtSciencePark Colours Scheme
###
### Website to convert HEX to RGB: http://hex.colorrrs.com.
### For some functions you should divide these numbers by 255.
###
###	No.	Color				HEX		RGB							CMYK					CHR		MAF/INFO
### --------------------------------------------------------------------------------------------------------------------
###	1	yellow				#FBB820 (251,184,32)				(0,26.69,87.25,1.57) 	=>	1 		or 1.0 > INFO
###	2	gold				#F59D10 (245,157,16)				(0,35.92,93.47,3.92) 	=>	2		
###	3	salmon				#E55738 (229,87,56) 				(0,62.01,75.55,10.2) 	=>	3 		or 0.05 < MAF < 0.2 or 0.4 < INFO < 0.6
###	4	darkpink			#DB003F ((219,0,63)					(0,100,71.23,14.12) 	=>	4		
###	5	lightpink			#E35493 (227,84,147)				(0,63,35.24,10.98) 		=>	5 		or 0.8 < INFO < 1.0
###	6	pink				#D5267B (213,38,123)				(0,82.16,42.25,16.47) 	=>	6		
###	7	hardpink			#CC0071 (204,0,113)					(0,0,0,0) 	=>	7		
###	8	lightpurple			#A8448A (168,68,138)				(0,0,0,0) 	=>	8		
###	9	purple				#9A3480 (154,52,128)				(0,0,0,0) 	=>	9		
###	10	lavendel			#8D5B9A (141,91,154)				(0,0,0,0) 	=>	10		
###	11	bluepurple			#705296 (112,82,150)				(0,0,0,0) 	=>	11		
###	12	purpleblue			#686AA9 (104,106,169)				(0,0,0,0) 	=>	12		
###	13	lightpurpleblue		#6173AD (97,115,173/101,120,180)	(0,0,0,0) 	=>	13		
###	14	seablue				#4C81BF (76,129,191)				(0,0,0,0) 	=>	14		
###	15	skyblue				#2F8BC9 (47,139,201)				(0,0,0,0) 	=>	15		
###	16	azurblue			#1290D9 (18,144,217)				(0,0,0,0) 	=>	16		 or 0.01 < MAF < 0.05 or 0.2 < INFO < 0.4
###	17	lightazurblue		#1396D8 (19,150,216)				(0,0,0,0) 	=>	17		
###	18	greenblue			#15A6C1 (21,166,193)				(0,0,0,0) 	=>	18		
###	19	seaweedgreen		#5EB17F (94,177,127)				(0,0,0,0) 	=>	19		
###	20	yellowgreen			#86B833 (134,184,51)				(0,0,0,0) 	=>	20		
###	21	lightmossgreen		#C5D220 (197,210,32)				(0,0,0,0) 	=>	21		
###	22	mossgreen			#9FC228 (159,194,40)				(0,0,0,0) 	=>	22		or MAF > 0.20 or 0.6 < INFO < 0.8
###	23	lightgreen			#78B113 (120,177,19)				(0,0,0,0) 	=>	23/X
###	24	green				#49A01D (73,160,29)					(0,0,0,0) 	=>	24/Y
###	25	grey				#595A5C (89,90,92)					(0,0,0,0) 	=>	25/XY	or MAF < 0.01 or 0.0 < INFO < 0.2
###	26	lightgrey			#A2A3A4	(162,163,164)				(0,0,0,0) 	=> 	26/MT
### 
### ADDITIONAL COLORS
### 27	midgrey				#D7D8D7
### 28	very lightgrey		#ECECEC
### 29	white				#FFFFFF
### 30	black				#000000
### --------------------------------------------------------------------------------------------------------------------

### ANNOTATION INFORMATION | FOR DEBUGGING
### ATHERO-EXPRESS DATA
### Annotations file Illumina Methylation 450K BeadChip
###
### Col.No. What                        Description                                                                                               Type  Example
### 1	      IlmnID	                    Unique CpG locus identifier from the Illumina CG database	                                                chr	  cg00035864 "cg00050873" "cg00061679" "cg00063477" ...
### 2	      Name	                      Unique CpG locus identifier from the Illumina CG database	                                                chr	  cg00035864 "cg00050873" "cg00061679" "cg00063477" ...
### 3	      AddressA_ID	                Address of probe A		      	      	      	      	                                                    int	  31729416 32735311 28780415 16712347 19779393 29674443 30703409 69792329 17691457 36709370 ...
### 4	      AlleleA_ProbeSeq	          Sequence for probe A                                                                                      chr	  AAAACACTAACAATCTTATCCACATAAACCCTTAAATTTATCTCAAATTC "ACAAAAAAACAACACACAACTATAATAATTTTTAAAATAAATAAACCCCA" "AAAACATTAAAAAACTAATTCACTACTATTTAATTACTTTATTTTCCATC" "TATTCTTCCACACAAAATACTAAACRTATATTTACAAAAATACTTCCATC" ...
### 5	      AddressB_ID	                Address of probe  B	                                                                                      int	  NA 31717405 NA NA NA 38703326 36767301 46723459 NA NA ...
### 6	      AlleleB_ProbeSeq	          Sequence for probe B	                                                                                    chr	  "ACGAAAAAACAACGCACAACTATAATAATTTTTAAAATAAATAAACCCCG" "" "" ...
### 7	      Infinium_Design_Type	      Defines Assay type - Infinium I or Infinium II	                                                          chr	  II "I" "II" "II" ...
### 8	      Next_Base	                  Base added at SBE step - Infinium I assays only	                                                          chr	  "A" "" "" ...
### 9	      Color_Channel	              Color of the incorporated base  (Red or Green) - Infinium I assays only	                                  chr	  "Red" "" "" ...
### 10	    Forward_Sequence	          Sequence (in 5'-3' orientation) flanking query site	                                                      chr	  AATCCAAAGATGATGGAGGAGTGCCCGCTCATGATGTGAAGTACCTGCTCAGCTGGAAAC[CG]AATTTGAGATAAATTCAAGGGTCTATGTGGACAAGACTGCTAGTGTCTCTCTCTGGATTG "TATCTCTGTCTGGCGAGGAGGCAACGCACAACTGTGGTGGTTTTTGGAGTGGGTGGACCC[CG]GCCAAGACGGCCTGGGCTGACCAGAGACGGGAGGCAGAAAAAGTGGGCAGGTGGTTGCAG" "TCAACAAATGAGAGACATTGAAGAACTAATTCACTACTATTTGGTTACTTTATTTTCCAT[CG]AAGAAAACCTCTTTTTAAAAACTAACACATAAATAAAATGAACGAAGAACAAACTAAACG" "CTCCTGTACTTGTTCATTAAATAATGATTCCTTGGATATACCAAGTCTGGATAGCGGATT[CG]ATGGAAGCATTTTTGTAAATATACGTTCAGTATTTTGTGTGGAAGAACACAATCTAGCTG" ...
### 11	    Genome_Build	              Genome build on which forward sequence is based	                                                          int	  37 37 37 37 37 37 37 37 37 37 ...
### 12	    CHR	                        Chromosome - genome build 37	                                                                            chr	  Y "Y" "Y" "Y" ...
### 13	    MAPINFO	                    Coordinates - genome build 37	                                                                            int	  8553009 9363356 25314171 22741795 21664296 21239348 8148233 15815688 4868996 6133740 ...
### 14	    SourceSeq	                  Unconverted design sequence	                                                                              chr	  AGACACTAGCAGTCTTGTCCACATAGACCCTTGAATTTATCTCAAATTCG "CGGGGTCCACCCACTCCAAAAACCACCACAGTTGTGCGTTGCCTCCTCGC" "CGATGGAAAATAAAGTAACCAAATAGTAGTGAATTAGTTCTTCAATGTCT" "CGATGGAAGCATTTTTGTAAATATACGTTCAGTATTTTGTGTGGAAGAAC" ...
### 15	    Chromosome_36	              Chromosome - genome build 36	                                                                            chr	  Y "Y" "Y" "Y" ...
### 16	    Coordinate_36	              Coordinates - genome build 36	                                                                            chr	  8613009 "9973356" "23723559" "21151183" ...
### 17	    Strand	                    Design strand	                                                                                            chr	  F "R" "R" "F" ...
### 18	    Probe_SNPs	                Assays with SNPs present within probe >10bp from query site	                                              chr	  "" "" "rs9341313" ...
### 19	    Probe_SNPs_10	              Assays with SNPs present within probe ≤10bp from query site (HM27 carryover or recently discovered)	      chr	  "" "" "rs13447379" ...
### 20	    Random_Loci	                Loci which were chosen randomly in the design proccess	                                                  logi  NA NA NA NA NA NA ...
### 21	    Methyl27_Loci	              Present or absent on HumanMethylation27 array		                                                          logi  NA NA NA NA NA NA ...
### 22	    UCSC_RefGene_Name	          Gene name (UCSC)		                      	                      	                                      chr	  TTTY18 "TSPY4;FAM197Y2" "DAZ1;DAZ4;DAZ4" "EIF1AY" ...
### 23	    UCSC_RefGene_Accession	    Accession number (UCSC)		                      	                      	                                chr	  NR_001550 "NM_001164471;NR_001553" "NM_004081;NM_020420;NM_001005375" "NM_004681" ...
### 24	    UCSC_RefGene_Group	        Gene region feature category (UCSC)		                      	                                            chr	  TSS1500 "Body;TSS1500" "Body;Body;Body" "Body" ...
### 25	    UCSC_CpG_Islands_Name	      CpG island name (UCSC)		                      	                      	                                chr	  "chrY:9363680-9363943" "" "chrY:22737825-22738052" ...
### 26	    Relation_to_UCSC_CpG_Island	Relationship to Canonical CpG Island: Shores - 0-2 kb from CpG island; Shelves - 2-4 kb from CpG island.	chr	  "N_Shore" "" "S_Shelf" ...
### 27	    Phantom	                    FANTOM-derived promoter	                                                                                  chr	  "" "" "" ...
### 28	    DMR	                        Differentially methylated region (experimentally determined)	                                            chr	  "" "" "" ...
### 29	    Enhancer	                  Enhancer element (informatically-determined)		      	                                                  logi  NA NA NA NA NA NA ...
### 30	    HMM_Island	                Hidden Markov Model Island		      	      	                                                            chr	  "Y:9973136-9976273" "" "" ...
### 31	    Regulatory_Feature_Name	    Regulatory feature (informatically determined)	                                                          chr	  "" "" "" ...
### 32	    Regulatory_Feature_Group	  Regulatory feature category		      	      	                                                            chr	  "" "" "" ...
### 33	    DHS	                        DNAse hypersensitive site (experimentally determined)	                                                    logi  NA NA NA NA NA NA ...
### 34	    UCSC_RefGene_Dist		        UCSC exact base pair position relative to UCSC Gene                                                       chr	  -1091 "127281;-480" "31067;-1665795;-1665795" "4199" ...
###
###
### After merging/cbinding RESULTS, STATS and ANNOTATIONS
### COLUMNS NOMINAL RESULTS EXAMPLE
### NO.	COLUMN-NAME             TYPE  Example
### 1.  ProbeID                       chr  "cg23024343" "cg23024343" "cg23024343" "cg23024343" ...
### 2.  VARIANT                       chr  "rs62483725" "rs2237678" "rs3779495" "rs2807" ...
### 3.  Distance_VARIANT_ProbeID      int  -176 1514 48509 59027 3918 15317 -57885 -53207 -53051 -39236 ...
### 4.  Nominal_P                     num  4.97e-57 4.97e-57 4.43e-57 1.45e-57 1.38e-56 ...
### 5.  Beta                          num  -0.087 -0.087 0.0865 0.0867 -0.0862 ...
### 6.  Z                             num  -15.9 -15.9 -15.9 -15.9 -15.8 ...
### 7.  SD                            num  0.00548 0.00548 -0.00545 -0.00543 0.00546 ...
### 8.  SEM                           num  0.00548 0.00548 -0.00545 -0.00544 0.00546 ...
### 9.  Bonferroni                    num  8.90e-51 8.90e-51 7.94e-51 2.60e-51 2.47e-50 ...
### 10. BenjHoch                     num  2.23e-51 2.23e-51 2.23e-51 2.23e-51 4.95e-51 ...
### 11. Q                            num  2.21e-51 2.21e-51 2.21e-51 2.21e-51 4.91e-51 ...
### 12. IlmnID                       chr  "cg23024343" "cg23024343" "cg23024343" "cg23024343" ...
### 13. ProbeID                      chr  "cg23024343" "cg23024343" "cg23024343" "cg23024343" ...
### 14. AddressA_ID                  int  47602338 47602338 47602338 47602338 47602338 47602338 47602338 47602338 47602338 47602338 ...
### 15. AlleleA_ProbeSeq             chr  "CATTCCTAAAAAATTAAACATTTCCAACAAAAACTTTAATTCACTATATC" "CATTCCTAAAAAATTAAACATTTCCAACAAAAACTTTAATTCACTATATC" "CATTCCTAAAAAATTAAACATTTCCAACAAAAACTTTAATTCACTATATC" "CATTCCTAAAAAATTAAACATTTCCAACAAAAACTTTAATTCACTATATC" ...
### 16. AddressB_ID                  int  NA NA NA NA NA NA NA NA NA NA ...
### 17. AlleleB_ProbeSeq             chr  "" "" "" "" ...
### 18. Infinium_Design_Type         chr  "II" "II" "II" "II" ...
### 19. Next_Base                    chr  "" "" "" "" ...
### 20. Color_Channel                chr  "" "" "" "" ...
### 21. Forward_Sequence             chr  "TACTATGCATCTGTGCTCACCAAGCCTATCAATCCGATACTCTGTCATTGGCCAAATCCC[CG]ATATAGTGAATCAAAGTTTCTGCTGGAAATGCTTAATTCCTCAGGAATGTTGGGGGTGGA" "TACTATGCATCTGTGCTCACCAAGCCTATCAATCCGATACTCTGTCATTGGCCAAATCCC[CG]ATATAGTGAATCAAAGTTTCTGCTGGAAATGCTTAATTCCTCAGGAATGTTGGGGGTGGA" "TACTATGCATCTGTGCTCACCAAGCCTATCAATCCGATACTCTGTCATTGGCCAAATCCC[CG]ATATAGTGAATCAAAGTTTCTGCTGGAAATGCTTAATTCCTCAGGAATGTTGGGGGTGGA" "TACTATGCATCTGTGCTCACCAAGCCTATCAATCCGATACTCTGTCATTGGCCAAATCCC[CG]ATATAGTGAATCAAAGTTTCTGCTGGAAATGCTTAATTCCTCAGGAATGTTGGGGGTGGA" ...
### 22. Genome_Build                 int  37 37 37 37 37 37 37 37 37 37 ...
### 23. CHR                          chr  "7" "7" "7" "7" ...
### 24. MAPINFO                      int  107201750 107201750 107201750 107201750 107201750 107201750 107201750 107201750 107201750 107201750 ...
### 25. SourceSeq                    chr  "ATTCCTGAGGAATTAAGCATTTCCAGCAGAAACTTTGATTCACTATATCG" "ATTCCTGAGGAATTAAGCATTTCCAGCAGAAACTTTGATTCACTATATCG" "ATTCCTGAGGAATTAAGCATTTCCAGCAGAAACTTTGATTCACTATATCG" "ATTCCTGAGGAATTAAGCATTTCCAGCAGAAACTTTGATTCACTATATCG" ...
### 26. Chromosome_36                chr  "7" "7" "7" "7" ...
### 27. Coordinate_36                chr  "106988986" "106988986" "106988986" "106988986" ...
### 28. Strand                       chr  "F" "F" "F" "F" ...
### 29. Probe_SNPs                   chr  "" "" "" "" ...
### 30. Probe_SNPs_10                chr  "" "" "" "" ...
### 31. Random_Loci                  logi  NA NA NA NA NA NA ...
### 32. Methyl27_Loci                logi  NA NA NA NA NA NA ...
### 33. UCSC_RefGene_Name            chr  "COG5;COG5;COG5" "COG5;COG5;COG5" "COG5;COG5;COG5" "COG5;COG5;COG5" ...
### 34. UCSC_RefGene_Accession       chr  "NM_006348;NM_001161520;NM_181733" "NM_006348;NM_001161520;NM_181733" "NM_006348;NM_001161520;NM_181733" "NM_006348;NM_001161520;NM_181733" ...
### 35. UCSC_RefGene_Group           chr  "Body;Body;Body" "Body;Body;Body" "Body;Body;Body" "Body;Body;Body" ...
### 36. UCSC_CpG_Islands_Name        chr  "chr7:107204114-107204797" "chr7:107204114-107204797" "chr7:107204114-107204797" "chr7:107204114-107204797" ...
### 37. Relation_to_UCSC_CpG_Island  chr  "N_Shelf" "N_Shelf" "N_Shelf" "N_Shelf" ...
### 38. Phantom                      chr  "" "" "" "" ...
### 39. DMR                          chr  "" "" "" "" ...
### 40. Enhancer                     logi  NA NA NA NA NA NA ...
### 41. HMM_Island                   chr  "" "" "" "" ...
### 42. Regulatory_Feature_Name      chr  "7:107201304-107201878" "7:107201304-107201878" "7:107201304-107201878" "7:107201304-107201878" ...
### 43. Regulatory_Feature_Group     chr  "Promoter_Associated_Cell_type_specific" "Promoter_Associated_Cell_type_specific" "Promoter_Associated_Cell_type_specific" "Promoter_Associated_Cell_type_specific" ...
### 44. DHS                          logi  NA NA NA NA NA NA ...
### 45. UCSC_RefGene_Dist            chr  "3208;3208;3208" "3208;3208;3208" "3208;3208;3208" "3208;3208;3208" ...
### 46. VARIANT                      chr  "rs62483725" "rs2237678" "rs3779495" "rs2807" ...
### 47. Chr                          int  7 7 7 7 7 7 7 7 7 7 ...
### 48. BP                           int  107201575 107203265 107250260 107260778 107205669 107217068 107143866 107148544 107148700 107162515 ...
### 49. OtherAlleleA                 chr  "G" "T" "G" "G" ...
### 50. CodedAlleleA                 chr  "A" "G" "A" "A" ...
### 51. MAF                          num  0.21 0.21 0.212 0.212 0.213 ...
### 52. MAC                          num  642 641 647 648 651 ...
### 53. CAF                          num  0.21 0.21 0.788 0.788 0.213 ...
### 54. AvgMAxPostCall               num  0.998 0.998 0.999 0.999 0.997 ...
### 55. Info                         num  0.995 0.995 0.998 0.997 0.994 ...
### 56. HWE                          num  0.2463 0.1885 0.1249 0.1459 0.0932 ...
### 57. N                            int  1526 1526 1526 1526 1526 1526 1526 1526 1526 1526 ...
### 58. Imputation                   chr  "imputed" "imputed" "imputed" "imputed" ...
###
### COLUMNS PERMUTATION RESULTS EXAMPLE
### NO.	COLUMN-NAME             TYPE  Example
### 1.  ProbeID                    : chr  "cg14590325" "cg23024343" "cg27284331" "cg19486437" ...
### 2.  NVariants                  : int  5442 5442 5442 5442 5442 5442 5442 5442 5442 5442 ...
### 3.  MLE_Beta_shape1            : num  1.039 0.911 1.033 1.002 1.041 ...
### 4.  MLE_Beta_shape2            : num  367 151 347 284 373 ...
### 5.  Dummy                      : num  369 284 363 343 373 ...
### 6.  VARIANT                    : chr  "rs2269778" "rs2807" "rs12705390" "rs2214916" ...
### 7.  Distance_VARIANT_ProbeID   : int  23153 59027 113087 -6764 18 3254 0 9317 259 -4201 ...
### 8.  Nominal_P                  : num  2.43e-44 1.45e-57 1.73e-36 5.56e-33 2.12e-28 ...
### 9.  Beta                       : num  -0.0444 0.0867 -0.0757 -0.0338 0.0571 ...
### 10.  Perm_P                     : num  1e-06 1e-06 1e-06 1e-06 1e-06 ...
### 11.  Approx_Perm_P              : num  2.08e-37 1.66e-33 1.17e-29 2.00e-24 2.49e-23 ...
### 12.  Z                          : num  -13.9 -15.9 -12.6 -11.9 -11 ...
### 13.  SD                         : num  0.00287 -0.00572 0.00567 0.00246 -0.0056 ...
### 14.  SEM                        : num  0.00319 -0.00544 0.00603 0.00284 -0.00519 ...
### 15.  Bonferroni                 : num  6.84e-35 5.47e-31 3.84e-27 6.59e-22 8.18e-21 ...
### 16.  BenjHoch                   : num  6.84e-35 2.74e-31 1.28e-27 1.65e-22 1.64e-21 ...
### 17.  Q                          : num  5.76e-35 2.30e-31 1.08e-27 1.39e-22 1.38e-21 ...
### 18.  IlmnID                     : chr  "cg14590325" "cg23024343" "cg27284331" "cg19486437" ...
### 19.  ProbeID                    : chr  "cg14590325" "cg23024343" "cg27284331" "cg19486437" ...
### 20.  AddressA_ID                : int  64795430 47602338 40603406 63657345 29749356 63760476 38669416 29707395 44666470 35736385 ...
### 21.  AlleleA_ProbeSeq           : chr  "AAATCTTAACTTCCACAAAAACCATTCCTAACAACCCTCTCTAAAATTTC" "CATTCCTAAAAAATTAAACATTTCCAACAAAAACTTTAATTCACTATATC" "AACRACACTTTACACTATCAACAATAAAACTTACCCTCAACAATTATAAC" "CTAACACTTAAATACTTACTAAATACCAAAAACTAAACTCAACAAATTAC" ...
### 22.  AddressB_ID                : int  NA NA NA NA NA NA NA 62623449 NA NA ...
### 23.  AlleleB_ProbeSeq           : chr  "" "" "" "" ...
### 24.  Infinium_Design_Type       : chr  "II" "II" "II" "II" ...
### 25.  Next_Base                  : chr  "" "" "" "" ...
### 26.  Color_Channel              : chr  "" "" "" "" ...
### 27.  Forward_Sequence           : chr  "TTTAATGTTTGCGGATTATTTTGAAGAAAATAAACGAGTGCTGTGGAGAGATTCATAGGA[CG]AAATTTCAGAGAGGGTTGTTAGGAATGGTCTTTGTGGAAGCTAAGATCTGATGGATGAGA" "TACTATGCATCTGTGCTCACCAAGCCTATCAATCCGATACTCTGTCATTGGCCAAATCCC[CG]ATATAGTGAATCAAAGTTTCTGCTGGAAATGCTTAATTCCTCAGGAATGTTGGGGGTGGA" "CACATCACAACAGCGACACTTTGCACTATCAACAATGAAGCTTGCCCTCAACAATTATGA[CG]TTACTGGTTTTAGTAACATAAAAATACATTGCTGGTTGACAGGAAGGATAAAAATGACAT" "CTTTTTCTAACAAAGGGATAATACTGCTATTCACCTTGTAGGCTTATTGTGACAATTTAA[CG]TAACCTGCTGAGCTCAGTTCCTGGTACTCAGCAAGCACTCAAGTGTTAGCTATGTTCTAG" ...
### 28.  Genome_Build               : int  37 37 37 37 37 37 37 37 37 37 ...
### 29.  CHR                        : chr  "7" "7" "7" "7" ...
### 30.  MAPINFO                    : int  107385716 107201750 106297689 108108892 106606094 107745446 107965081 106694832 108097477 106592626 ...
### 31.  SourceSeq                  : chr  "GATCTTAGCTTCCACAAAGACCATTCCTAACAACCCTCTCTGAAATTTCG" "ATTCCTGAGGAATTAAGCATTTCCAGCAGAAACTTTGATTCACTATATCG" "CGTCATAATTGTTGAGGGCAAGCTTCATTGTTGATAGTGCAAAGTGTCGC" "CGTAACCTGCTGAGCTCAGTTCCTGGTACTCAGCAAGCACTCAAGTGTTA" ...
### 32.  Chromosome_36              : chr  "7" "7" "7" "7" ...
### 33.  Coordinate_36              : chr  "107172952" "106988986" "106084925" "107896128" ...
### 34.  Strand                     : chr  "F" "F" "R" "F" ...
### 35.  Probe_SNPs                 : chr  "" "" "" "" ...
### 36.  Probe_SNPs_10              : chr  "" "" "" "" ...
### 37.  Random_Loci                : logi  NA NA NA NA NA TRUE ...
### 38.  Methyl27_Loci              : logi  NA NA NA NA NA NA ...
### 39.  UCSC_RefGene_Name          : chr  "CBLL1;CBLL1" "COG5;COG5;COG5" "" "" ...
### 40.  UCSC_RefGene_Accession     : chr  "NM_024814;NR_024199" "NM_006348;NM_001161520;NM_181733" "" "" ...
### 41.  UCSC_RefGene_Group         : chr  "Body;Body" "Body;Body;Body" "" "" ...
### 42.  UCSC_CpG_Islands_Name      : chr  "chr7:107383657-107385021" "chr7:107204114-107204797" "chr7:106300402-106301573" "" ...
### 43.  Relation_to_UCSC_CpG_Island: chr  "S_Shore" "N_Shelf" "N_Shelf" "" ...
### 44.  Phantom                    : chr  "" "" "" "" ...
### 45.  DMR                        : chr  "RDMR" "" "" "" ...
### 46.  Enhancer                   : logi  NA NA NA TRUE TRUE NA ...
### 47.  HMM_Island                 : chr  "" "" "" "" ...
### 48.  Regulatory_Feature_Name    : chr  "" "7:107201304-107201878" "" "" ...
### 49.  Regulatory_Feature_Group   : chr  "" "Promoter_Associated_Cell_type_specific" "" "" ...
### 50.  DHS                        : logi  NA NA NA NA NA NA ...
### 51.  UCSC_RefGene_Dist          : chr  "1575;1575" "3208;3208;3208" "" "" ...
### 52.  VARIANT                    : chr  "rs2269778" "rs2807" "rs12705390" "rs2214916" ...
### 53.  Chr                        : int  7 7 7 7 7 7 7 7 7 7 ...
### 54.  BP                         : int  107408870 107260778 106410777 108102129 106606113 107748701 107965082 106704150 108097737 106588426 ...
### 55.  OtherAlleleA               : chr  "G" "G" "G" "A" ...
### 56.  CodedAlleleA               : chr  "A" "A" "A" "T" ...
### 57.  MAF                        : num  0.448 0.212 0.218 0.383 0.237 ...
### 58.  MAC                        : num  1368 648 667 1170 725 ...
### 59.  CAF                        : num  0.448 0.788 0.218 0.383 0.237 ...
### 60.  AvgMAxPostCall             : num  0.997 0.999 0.999 0.997 0.991 ...
### 61.  Info                       : num  0.995 0.997 0.998 0.996 0.982 ...
### 62.  HWE                        : num  0.255 0.146 0.881 0.745 0.671 ...
### 63.  N                          : int  1526 1526 1526 1526 1526 1526 1526 1526 1526 1526 ...
### 64.  Imputation                 : chr  "imputed" "imputed" "genotyped" "imputed" ...
###
### CTMMM DATA
###
### After merging/cbinding RESULTS, STATS and ANNOTATIONS
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
### NO.	COLUMN-NAME             TYPE  Example
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


