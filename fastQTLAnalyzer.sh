#!/bin/bash
#
# You can use the variables below (indicated by "#$") to set some things for the 
# submission system.
#$ -S /bin/bash # the type of BASH you'd like to use
#$ -N fastQTLAnalyzer_v2 # the name of this script
#$ -hold_jid some_other_basic_bash_script # the current script (basic_bash_script) will hold until some_other_basic_bash_script has finished
#$ -o /hpc/dhl_ec/svanderlaan/projects/test_mqtl/fastQTLAnalyzer_v2.log # the log file of this job
#$ -e /hpc/dhl_ec/svanderlaan/projects/test_mqtl/fastQTLAnalyzer_v2.errors # the error file of this job
#$ -l h_rt=04:00:00 # h_rt=[max time, hh:mm:ss, e.g. 02:02:01] - this is the time you think the script will take
#$ -l h_vmem=8G #  h_vmem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
# -l tmpspace=32G # this is the amount of temporary space you think your script will use
#$ -M s.w.vanderlaan-2@umcutrecht.nl # you can send yourself emails when the job is done; "-M" and "-m" go hand in hand
#$ -m ea # you can choose: b=begin of job; e=end of job; a=abort of job; s=suspended job; n=no mail is send
#$ -cwd # set the job start to the current directory - so all the things in this script are relative to the current directory!!!
#
# Another useful tip: you can set a job to run after another has finished. Name the job 
# with "-N SOMENAME" and hold the other job with -hold_jid SOMENAME". 
# Further instructions: https://wiki.bioinformatics.umcutrecht.nl/bin/view/HPC/HowToS#Run_a_job_after_your_other_jobs
#
# It is good practice to properly name and annotate your script for future reference for
# yourself and others. Trust me, you'll forget why and how you made this!!!

### REGARDING NOTES ###
### Please note that uncommented notes can be found at the end of this script.
###

### MoSCoW FEATURE LIST ###
###
### * source arguments list, so ./fastQTLAnalyzer.sh arguments_file.txt
###
### --- THESE COULD BE ARGUMENTS --- ###
### * REGIONS=${PROJECT}/regions_for_qtl.small.txt # regions_for_eqtl.txt OR regions_for_qtl.small.txt; [arg5]
### * PROJECTNAME="CAD" # [arg6]
### * EXCLUSION_COV="${PROJECT}/excl_cov.txt" # [arg7]
### * EXCLUSION_TYPE="DEFAULT" # DEFAULT/SMOKER/NONSMOKER/MALES/FEMALES/T2D/NONT2D [CKD/NONCKD/PRE2007/POST2007/NONAEGS/NONAEGSFEMALES/NONAEGSMALES -- these are AE-specific!!!]  # [arg8]
### * 
### * EMAIL="s.w.vanderlaan-2@umcutrecht.nl" # [arg9]
### * MAILTYPE="as" # [arg10]
### * 
### * PERMSTART="100" # [arg9]
### * PERMEND="1000" # [arg7]
### * MAF="0.005" # [arg8]
### * INFO="0.9" # [arg9]
### * HWE="6" # [arg10]
### * 
### * QUEUE_QCTOOL="01:00:00" # [arg11]
### * VMEM_QCTOOL="4G" # [arg12]
### * QUEUE_NOM="08:00:00" # [arg13]
### * VMEM_NOM="4G" # [arg14]
### * QUEUE_PERM="24:00:00" # [arg15]
### * VMEM_PERM="4G" # [arg16]
### * 
### * 
### * 
### * 

script_copyright_message() {
	echo ""
	echo ""
	echo ""
	THISYEAR=$(date +'%Y')
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "+ The MIT License (MIT)                                                                                 +"
	echo "+ Copyright (c) 2015-${THISYEAR} Sander W. van der Laan                                                        +"
	echo "+                                                                                                       +"
	echo "+ Permission is hereby granted, free of charge, to any person obtaining a copy of this software and     +"
	echo "+ associated documentation files (the \"Software\"), to deal in the Software without restriction,         +"
	echo "+ including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, +"
	echo "+ and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, +"
	echo "+ subject to the following conditions:                                                                  +"
	echo "+                                                                                                       +"
	echo "+ The above copyright notice and this permission notice shall be included in all copies or substantial  +"
	echo "+ portions of the Software.                                                                             +"
	echo "+                                                                                                       +"
	echo "+ THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT     +"
	echo "+ NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND                +"
	echo "+ NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES  +"
	echo "+ OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN   +"
	echo "+ CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                            +"
	echo "+                                                                                                       +"
	echo "+ Reference: http://opensource.org.                                                                     +"
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
}

script_arguments_error() {
	echo "$1" # additional error message
	echo ""
	echo " * Argument #1  indicate which study type you want to analyze, so either:"
	echo "                [AEMS450K1/AEMS450K2/CTMM]:"
	echo "                - AEMS450K1: methylation quantitative trait locus (mQTL) analysis "
	echo "                             on plaques or blood in the Athero-Express Methylation "
	echo "                             Study 450K phase 1."
	echo "                - AEMS450K2: mQTL analysis on plaques or blood in the Athero-Express"
	echo "                             Methylation Study 450K phase 2."
	echo "                - CTMM:      expression QTL (eQTL) analysis in monocytes from CTMM."
	echo ""
	echo " * Argument #2  the sample type must be [AEMS450K1: PLAQUES/BLOOD], "
	echo "                [AEMS450K2: PLAQUES], or [CTMM: MONOCYTES]."
	echo ""
	echo " * Argument #3  the root directory, e.g. /hpc/dhl_ec/svanderlaan/projects/test_qtl."
	echo " "
	echo " * Argument #4  where you want stuff to be save inside the rootdir, "
	echo "                e.g. mqtl_aems450k1."
	echo ""
	echo " An example command would be: "
	echo ""
	echo "./fastQTLAnalyzer.sh [arg1] [arg2] [arg3] [arg4]"
	echo ""
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  	# The wrong arguments are passed, so we'll exit the script now!
  	script_copyright_message
  	exit 1
}

echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "+                                      QUANTITATIVE TRAIT LOCUS ANALYZER                                +"
echo "+                                                                                                       +"
echo "+                                                                                                       +"
echo "+ * Written by  : Sander W. van der Laan                                                                +"
echo "+ * E-mail      : s.w.vanderlaan-2@umcutrecht.nl                                                        +"
echo "+ * Last update : 2016-10-28                                                                            +"
echo "+ * Version     : 2.0.5                                                                                 +"
echo "+                                                                                                       +"
echo "+ * Description : This script will set some directories, execute something in a for-loop, and will then +"
echo "+                 submit this in a job.                                                                 +"
echo "+                                                                                                       +"
echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Today's: "$(date)
TODAY=$(date +"%Y%m%d")
echo ""

### SET STUDY AND SAMPLE TYPE
### Note: All analyses with AE data are presumed to be constrained to CEA-patients only.
###       You can set the exclusion criteria 'NONAEGS/FEMALES/MALES' if you want to analyse
###       all AE data!
### Set the analysis type.
STUDY_TYPE=${1} # AEMS450K1/AEMS450K2/CTMM

### Set the analysis type.
SAMPLE_TYPE=${2} # AE: PLAQUES/BLOOD; CTMM: MONOCYTES

### START of if-else statement for the number of command-line arguments passed ###
if [[ $# -lt 4 ]]; then 
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "                        *** Oh no! Computer says no! ***"
	echo ""
	echo "Number of arguments found "$#"."
	echo " [arg1]: ${STUDY_TYPE}"
	echo " [arg2]: ${SAMPLE_TYPE}"
	echo " [arg3]: ${3}" #ROOTDIR
	echo " [arg4]: ${4}" #PROJECTDIR with ROOTDIR
	echo ""
	script_arguments_error "You must supply at least [4] arguments when running a mQTL or eQTL analysis using Athero-Express or CTMM data!"
	
elif [[ (${STUDY_TYPE} = "AEMS450K1" || ${STUDY_TYPE} = "AEMS450K2") && ${SAMPLE_TYPE} = "MONOCYTES" ]]; then 
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "                        *** Oh no! Computer says no! ***"
	echo ""
	echo "Number of arguments found "$#"."
	echo " [arg1]: ${STUDY_TYPE}"
	echo " [arg2]: ${SAMPLE_TYPE}"
	echo " [arg3]: ${3}" #ROOTDIR
	echo " [arg4]: ${4}" #PROJECTDIR with ROOTDIR
	echo ""
	script_arguments_error "When running a *** mQTL analysis *** using ${STUDY_TYPE}, you must supply 'PLAQUES' or 'BLOOD' as SAMPLE_TYPE [arg2]!"
	
elif [[ ${STUDY_TYPE} = "AEMS450K2" && ${SAMPLE_TYPE} = "BLOOD" ]]; then 
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "                        *** Oh no! Computer says no! ***"
	echo ""
	echo "Number of arguments found "$#"."
	echo " [arg1]: ${STUDY_TYPE}"
	echo " [arg2]: ${SAMPLE_TYPE}"
	echo " [arg3]: ${3}" #ROOTDIR
	echo " [arg4]: ${4}" #PROJECTDIR with ROOTDIR
	echo ""
	script_arguments_error "When running a *** mQTL analysis *** using ${STUDY_TYPE}, you must supply 'PLAQUES' as SAMPLE_TYPE [arg2]!"
	
elif [[ ${STUDY_TYPE} = "CTMM" && (${SAMPLE_TYPE} = "PLAQUES" || ${SAMPLE_TYPE} = "BLOOD") ]]; then 
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "                        *** Oh no! Computer says no! ***"
	echo ""
	echo "Number of arguments found "$#"."
	echo " [arg1]: ${STUDY_TYPE}"
	echo " [arg2]: ${SAMPLE_TYPE}"
	echo " [arg3]: ${3}" #ROOTDIR
	echo " [arg4]: ${4}" #PROJECTDIR with ROOTDIR
	echo ""
	script_arguments_error "When running a *** eQTL analysis *** using ${STUDY_TYPE}, you must supply 'MONOCYTES' as SAMPLE_TYPE [arg2]!"
	
else
	
	### GENERIC SETTINGS
	SOFTWARE=/hpc/local/CentOS7/dhl_ec/software
	QCTOOL=${SOFTWARE}/qctool_v1.5-linux-x86_64-static/qctool
	SNPTEST252=${SOFTWARE}/snptest_v2.5.2_CentOS6.5_x86_64_static/snptest_v2.5.2
	FASTQTL=${SOFTWARE}/fastqtl_v2.184
	FASTQCTLADDON=${SOFTWARE}/fastqtl_supportives
	FASTQTLPARSER=${FASTQCTLADDON}/NominalResultsParser.py
	LZ13=${SOFTWARE}/locuszoom_1.3/bin/locuszoom
	BGZIP=${SOFTWARE}/htslib-1.3/bgzip
	TABIX=${SOFTWARE}/htslib-1.3/tabix

	### PROJECT SPECIFIC -- PLEASE CHANGE TO YOUR SITUATION
	### --- THESE COULD BE ARGUMENTS --- ###
	ROOTDIR=${3} # the root directory, e.g. /hpc/dhl_ec/svanderlaan/projects/test_qtl; [arg3]
	PROJECTDIR=${4} # where you want stuff to be save inside the rootdir, e.g. mqtl_aems450k1; [arg4]
	if [ ! -d ${ROOTDIR}/${PROJECTDIR} ]; then
				echo "The project directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${ROOTDIR}/${PROJECTDIR}
			else
				echo "The project directory '${PROJECTDIR}' already exists."
			fi
	PROJECT=${ROOTDIR}/${PROJECTDIR}
	
	### FASTQTL & QCTOOL SETTINGS
	SEEDNO=91149216
	### --- THESE COULD BE ARGUMENTS --- ###
	REGIONS=${ROOTDIR}/regions_for_qtl.small.txt # regions_for_eqtl.txt OR regions_for_qtl.small.txt; [arg5]
	PERMSTART="100" # [arg6]
	PERMEND="1000" # [arg7]
	MAF="0.005" # [arg8]
	INFO="0.9" # [arg9]
	HWE="6" # [arg10]

	### SETTING STUDY AND SAMPLE TYPE SPECIFIC THINGS
	if [[ ${STUDY_TYPE} == "AEMS450K1" ]]; then
		### AEMS450K1 SPECIFIC -- DO NOT CHANGE
		ORIGINALS=/hpc/dhl_ec/data/_ae_originals
		GENETICDATA=${ORIGINALS}/AEGS_COMBINED_IMPUTE2_1000Gp3_GoNL5
		AEMS450K1=${ORIGINALS}/AEMethylS_IlluminaMethylation450K
	
		### for file names
		STUDYNAME="aegs"
		STUDYJOBNAME="AEMS450K1"
		SNPTESTDATA="aegs_combo_1kGp3GoNL5_RAW_chr"
		#SNPTESTSAMPLEDATA=""
		SNPTESTOUTPUTDATA="aegs_1kGp3GoNL5"
		if [[ ${SAMPLE_TYPE} == "PLAQUES" ]]; then
			FASTQTLDATA="${AEMS450K1}/AEM_mQTL_INPUT_DATA/aems450k1_QC_443872_plaques.bed.gz"
			FASTQTLINDEX="${AEMS450K1}/AEM_mQTL_INPUT_DATA/aems450k1_QC_443872_plaques.bed.gz.tbi"
		elif [[ ${STUDY_TYPE} == "BLOOD" ]]; then
			FASTQTLDATA="${AEMS450K1}/AEM_mQTL_INPUT_DATA/aems450k1_QC_443872_blood.bed.gz"
			FASTQTLINDEX="${AEMS450K1}/AEM_mQTL_INPUT_DATA/aems450k1_QC_443872_blood.bed.gz.tbi"
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'sample type' please: '${SAMPLE_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
		### COVARIATES FILE
		COVARIATES="${GENETICDATA}/covariates_aegs_combo.all.cov"
	
	elif [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
		### AEMS450K2 SPECIFIC -- DO NOT CHANGE
		ORIGINALS=/hpc/dhl_ec/data/_ae_originals
		GENETICDATA=${ORIGINALS}/AEGS_COMBINED_IMPUTE2_1000Gp3_GoNL5
		AEMS450K2=${ORIGINALS}/AEMS450K2
	
		### for file names
		STUDYNAME="aegs"
		STUDYJOBNAME="AEMS450K2"
		SNPTESTDATA="aegs_combo_1kGp3GoNL5_RAW_chr"
		#SNPTESTSAMPLEDATA=""
		SNPTESTOUTPUTDATA="aegs_1kGp3GoNL5"
		if [[ ${SAMPLE_TYPE} == "PLAQUES" ]]; then
			FASTQTLDATA="${AEMS450K2}/AEM_mQTL_INPUT_DATA/aems450k2_QC_443872_plaques.bed.gz"
			FASTQTLINDEX="${AEMS450K2}/AEM_mQTL_INPUT_DATA/aems450k2_QC_443872_plaques.bed.gz.tbi"
		elif [[ ${STUDY_TYPE} == "BLOOD" ]]; then
			FASTQTLDATA="${AEMS450K2}/AEM_mQTL_INPUT_DATA/aems450k2_QC_443872_blood.bed.gz"
			FASTQTLINDEX="${AEMS450K2}/AEM_mQTL_INPUT_DATA/aems450k2_QC_443872_blood.bed.gz.tbi"
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'sample type' please: '${SAMPLE_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
		### COVARIATES FILE
		COVARIATES="${GENETICDATA}/covariates_aegs_combo.all.cov"
	
	elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
		### CTMM SPECIFIC -- DO NOT CHANGE
		ORIGINALS=/hpc/dhl_ec/data/_ctmm_originals
		GENETICDATA=${ORIGINALS}/CTMMAxiomTX_IMPUTE2_1000Gp3_GoNL5
		CTMMEXPRESSIONDATA=${ORIGINALS}/CTMMHumanHT12v4r2_15002873B
	
		### for file names
		STUDYNAME="ctmm"
		STUDYJOBNAME="CTMM"
		SNPTESTDATA="ctmm_1kGp3GoNL5_RAW_chr"
		#SNPTESTSAMPLEDATA=""
		SNPTESTOUTPUTDATA="ctmm_1kGp3GoNL5"
		if [[ ${SAMPLE_TYPE} == "MONOCYTES" ]]; then
			FASTQTLDATA="${CTMMEXPRESSIONDATA}/phenotype_ctmm.all.bed.gz"
			FASTQTLINDEX="${CTMMEXPRESSIONDATA}/phenotype_ctmm.all.bed.gz.tbi"
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'sample type' please: '${SAMPLE_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
		### COVARIATES FILE
		COVARIATES="${GENETICDATA}/covariates_ctmm.all.cov"
	else
		echo "                        *** ERROR *** "
		echo "Something is rotten in the City of Gotham; most likely a typo. "
		echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
		echo "                *** END OF ERROR MESSAGE *** "
		exit 1
	fi


	### THIS CODE IS PROBABLY NOT NEEDED -- IT IS INTENDED AS ANALYSIS TYPE FOR-LOOP 
	### AE-SPECIFIC -- see also notes under 'SET STUDY TYPE'
	###for TYPE in DEFAULT MALES FEMALES SMOKER NONSMOKER T2D NONT2D CKD NONCKD PRE2007 POST2007 NONAEGS NONAEGSFEMALES NONAEGSMALES; do
	### CTMM SPECIFIC
	###for TYPE in DEFAULT MALES FEMALES SMOKER NONSMOKER T2D NONT2D; do
	### CODE ABOVE PROBABLY OBSOLETE

	### SET EXCLUSION TYPE
	### --- THESE COULD BE ARGUMENTS --- ###
	EXCLUSION_TYPE="DEFAULT" # DEFAULT/SMOKER/NONSMOKER/MALES/FEMALES/T2D/NONT2D [CKD/NONCKD/PRE2007/POST2007/NONAEGS/NONAEGSFEMALES/NONAEGSMALES -- these are AE-specific!!!]

	### QSUB SETTINGS
	### --- THESE COULD BE ARGUMENTS --- ###
	QUEUE_QCTOOL="01:00:00"
	VMEM_QCTOOL="4G"
	QUEUE_NOM="08:00:00"
	VMEM_NOM="4G"
	QUEUE_PERM="24:00:00"
	VMEM_PERM="4G"

	### MAIL SETTINGS -- PLEASE CHANGE TO YOUR SITUATION
	### --- THESE COULD BE ARGUMENTS --- ###
	EMAIL="s.w.vanderlaan-2@umcutrecht.nl"
	MAILTYPE="as"

	### PROJECT SPECIFIC -- PLEASE CHANGE TO YOUR SITUATION
	### --- THESE COULD BE ARGUMENTS --- ###
	PROJECTNAME="CAD"
	EXCLUSION_COV="${ROOTDIR}/excl_cov.txt"

	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "The following is set:"
	echo ""
	echo "Software directory                                 ${SOFTWARE}"
	echo "Where \"qctool\" resides                             ${QCTOOL}"
	echo "Where \"fastQTL\" resides                            ${FASTQTL}"
	echo "Where \"bgzip\" resides                              ${BGZIP}"
	echo "Where \"tabix\" resides                              ${TABIX}"
	echo "Where \"snptest 2.5.2\" resides                      ${SNPTEST252}"
	echo ""

	echo "Original Athero-Express/CTMM data directory        ${ORIGINALS}"
	echo "AEGS/CTMM genetic data directory (1kGp3v5+GoNL5)   ${GENETICDATA}"
	echo ""
	echo "Expression or methylation data directory           ${CTMMEXPRESSIONDATA}${AEMS450K1}${AEMS450K2}"
	echo ""     
	echo "Project directory                                  ${PROJECT}"
	echo ""     
	echo "Additional fastQTL specific settings:"     
	echo ""     
	echo "Seed number                                        ${SEEDNO}"
	echo ""     
	echo "We will run this script on                         ${TODAY}"
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

	### AEMS - CTMM SPECIFIC
	### CREATES RESULTS AND SUMMARY DIRECTORIES; SETS EXCLUSION-FILES
	### DEFAULT
	if [[ ${EXCLUSION_TYPE} == "DEFAULT" ]]; then
		# Results directory
		if [ ! -d ${PROJECT}/qtl ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl
		# Summary directory
		if [ ! -d ${PROJECT}/qtl_summary ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_forFastQTL.list"

		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
	### MALES ONLY ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "FEMALES" ]]; then
		# Results directory
		if [ ! -d ${PROJECT}/qtl_males ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_males
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_males
		# Summary directory
		if [ ! -d ${PROJECT}/qtl_summary_males ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_males
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_males
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA_Females.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA_Females.list"

		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_FEMALES.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_FEMALES_forFastQTL.list"
			
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
	
	### FEMALES ONLY ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "MALES" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_females ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_females
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_females
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_females ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_females
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_females
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA_Males.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA_Males.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_MALES.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_MALES_forFastQTL.list"
			
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi

	### SMOKER ONLY ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "NONSMOKER" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_smoker ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_smoker
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_smoker
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_smoker ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_smoker
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_smoker
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA_nonSMOKER.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA_nonSMOKER.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_nonSMOKER.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_nonSMOKER_forFastQTL.list"
			
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
	
	### NONSMOKER ONLY ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "SMOKER" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_nonsmoker ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_nonsmoker
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_nonsmoker
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_nonsmoker ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_nonsmoker
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_nonsmoker
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA_SMOKER.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA_SMOKER.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_SMOKER.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_SMOKER_forFastQTL.list"
			
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
	
	### NON-TYPE 2 DIABETES ONLY ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "T2D" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_nont2d ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_nont2d
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_nont2d
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_nont2d ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_nont2d
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_nont2d
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA_T2D.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA_T2D.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_T2D.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_T2D_forFastQTL.list"
			
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
	
	### TYPE 2 DIABETES ONLY ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "NONT2D" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_t2d ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_t2d
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_t2d
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_t2d ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_t2d
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_t2d
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA_nonT2D.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA_nonT2D.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_nonT2D.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_nonT2D_forFastQTL.list"
			
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi

	##### THIS PART IS ATHERO-EXPRESS SPECIFIC ONLY #####

	### NON-CKD ONLY ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "CKD" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_nonckd ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_nonckd
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_nonckd
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_nonckd ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_nonckd
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_nonckd
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA_CKD.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA_CKD.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			echo "The exclusion criterium '${EXCLUSION_TYPE}' does *not* exist for CTMM."
			exit 1
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi

	### CKD DIABETES ONLY ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "NONCKD" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_ckd ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_ckd
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_ckd
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_ckd ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_ckd
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_ckd
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA_nonCKD.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA_nonCKD.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			echo "The exclusion criterium '${EXCLUSION_TYPE}' does *not* exist for CTMM."
			exit 1
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi

	### POST 2007 ONLY ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "PRE2007" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_post2007 ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_post2007
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_post2007
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_post2007 ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_post2007
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_post2007
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA_pre2007.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA_pre2007.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			echo "The exclusion criterium '${EXCLUSION_TYPE}' does *not* exist for CTMM."
			exit 1
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi

	### PRE 2007 ONLY ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "POST2007" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_pre2007 ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_pre2007
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_pre2007
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_pre2007 ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_pre2007
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_pre2007
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCEA_post2007.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCEA_post2007.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			echo "The exclusion criterium '${EXCLUSION_TYPE}' does *not* exist for CTMM."
			exit 1
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi

	### ALL AEGS ANALYSIS
	elif [[ ${EXCLUSION_TYPE} == "NONAEGS" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_allaegs ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_allaegs
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_allaegs
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_allaegs ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_allaegs
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_allaegs
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonAEGS.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonAEGS.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			echo "The exclusion criterium '${EXCLUSION_TYPE}' does *not* exist for CTMM."
			exit 1
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
	
	### ALL AEGS ANALYSIS -- MALES ONLY
	elif [[ ${EXCLUSION_TYPE} == "NONAEGSFEMALES" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_allaegs_males ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_allaegs_males
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_allaegs_males
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_allaegs_males ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_allaegs_males
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_allaegs_males
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonFemales.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonFemales.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			echo "The exclusion criterium '${EXCLUSION_TYPE}' does *not* exist for CTMM."
			exit 1
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
	
	### ALL AEGS ANALYSIS -- FEMALES ONLY
	elif [[ ${EXCLUSION_TYPE} == "NONAEGSMALES" ]]; then		
		# Results directory
		if [ ! -d ${PROJECT}/qtl_allaegs_females ]; then
				echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
				mkdir -v ${PROJECT}/qtl_allaegs_females
			else
				echo "The regional directory already exists."
			fi
		RESULTS=${PROJECT}/qtl_allaegs_females
		# Summary directory	
		if [ ! -d ${PROJECT}/qtl_summary_allaegs_females ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/qtl_summary_allaegs_females
		else
			echo "The regional directory already exists."
		fi
		SUMMARY=${PROJECT}/qtl_summary_allaegs_females
	
		# Exclusion files
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonMales.list"
			EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonMales.list"
			
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			echo "The exclusion criterium '${EXCLUSION_TYPE}' does *not* exist for CTMM."
			exit 1
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi

	else
		echo "                        *** ERROR *** "
		echo "Something is rotten in the City of Gotham; most likely a typo. "
		echo "Double back, and check you 'exclusion type' please."	
		echo "                *** END OF ERROR MESSAGE *** "
		exit 1
	fi

	### OVERVIEW OF REGIONS
	echo ""
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "The list of regions to investigate:"
		echo "* Variant	Locus	Chr	BP	StartRange	EndRange	WindowSize	Type	Phenotype"
		while IFS='' read -r REGION || [[ -n "$REGION" ]]; do
		LINE=${REGION}
		echo "* ${LINE}"
		done < ${REGIONS}
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

	### EXTRACTION OF DATA AND ANALYSIS
	echo ""
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "Extracting loci with the specified bp range."
	while IFS='' read -r REGIONOFINTEREST || [[ -n "$REGIONOFINTEREST" ]]; do
		###	1		2		3	4	5		6		7			8		9
		###Variant	Locus	Chr	BP	StartRange	EndRange	WindowSize	Type	Phenotype
		LINE=${REGIONOFINTEREST}
		VARIANT=$(echo "${LINE}" | awk '{print $1}')
		LOCUS=$(echo "${LINE}" | awk '{print $2}')
		CHR=$(echo "${LINE}" | awk '{print $3}')
		BP=$(echo "${LINE}" | awk '{print $4}')
		START=$(echo "${LINE}" | awk '{print $5}')
		END=$(echo "${LINE}" | awk '{print $6}')
		WINDOWSIZE=$(echo "${LINE}" | awk '{print $7}')
		TYPE=$(echo "${LINE}" | awk '{print $8}')
		PHENOTYPE=$(echo "${LINE}" | awk '{print $9}')
	
		echo ""
		echo ""
		echo "========================================================================================================="
		echo "Processing ${VARIANT} locus on ${CHR} between ${START} and ${END}..."
		echo "========================================================================================================="
		### Make directories for script if they do not exist yet (!!!PREREQUISITE!!!)
		if [ ! -d ${RESULTS}/${VARIANT}_${PROJECTNAME} ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${RESULTS}/${VARIANT}_${PROJECTNAME}
		else
			echo "The regional directory already exists."
		fi
		
		REGIONALDIR=${RESULTS}/${VARIANT}_${PROJECTNAME}
		
		### Extraction relevant regions for QTL analysis using fastQTL
		echo ""
		echo "* Creating bash-script to submit qctool extraction of region ${CHR}:${START}-${END} near ${LOCUS}..."
		### for DEBUGGING
		${QCTOOL} -g ${GENETICDATA}/${SNPTESTDATA}${CHR}.bgen -s ${GENETICDATA}/${SNPTESTDATA}${CHR}.sample -og ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -os ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -incl-range ${CHR}:${START}-${END}
		#echo "${QCTOOL} -g ${GENETICDATA}/${SNPTESTDATA}${CHR}.bgen -s ${GENETICDATA}/${SNPTESTDATA}${CHR}.sample -og ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -os ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -incl-range ${CHR}:${START}-${END} "> ${REGIONALDIR}/${STUDYNAME}_genex_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
		#qsub -S /bin/bash -N GENEX${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_genex_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_genex_${VARIANT}_excl_${EXCLUSION_TYPE}.log -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_genex_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
		
		### Applying some QC metrics on the extracted data -- exclude samples: -excl-samples ${EXCLUSION_NONAEMS450K1} 
		echo "* Creating bash-script to submit qctool filtering of the ${LOCUS} based on MAF > ${MAF}, INFO > ${INFO} and HWE -log10(p) > ${HWE}..."
		### for DEBUGGING
		${QCTOOL} -g ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -s ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -og ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -os ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -maf ${MAF} 1 -info ${INFO} 1 -hwe ${HWE} 
		#echo "${QCTOOL} -g ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -s ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -og ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -os ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -maf ${MAF} 1 -info ${INFO} 1 -hwe ${HWE} "> ${REGIONALDIR}/${STUDYNAME}_genqc_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
		#qsub -S /bin/bash -N GENQC${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid GENEX${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_genqc_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_genqc_${VARIANT}_excl_${EXCLUSION_TYPE}.log -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_genqc_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
		
		#### Calculating statistics 
		echo "* Creating bash-script to submit to calculate summary statistics of region ${CHR}:${START}-${END}..."
		### for DEBUGGING -- to exclude samples: -exclude_samples ${EXCLUSION_NONAEMS450K1} 
		${SNPTEST252} -data ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -summary_stats_only -hwe -o ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats 
		#echo "${SNPTEST252} -data ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -summary_stats_only -hwe -o ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats "> ${REGIONALDIR}/${STUDYNAME}_genstats_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
		#qsub -S /bin/bash -N GENSTAT${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid GENQC${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_genstats_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_genstats_${VARIANT}_excl_${EXCLUSION_TYPE}.log -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_genstats_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
		
		#### Make VCF
		#### example: qctool_v15 -g ctmm_1kGp3GoNL5_QC_chr7.7q22.gen.gz -s ctmm_phenocov.sample -og ctmm_1kGp3GoNL5_QC_chr7.7q22.vcf
		echo "* Creating bash-script to submit VCF-file generation..."
		### for DEBUGGING
		${QCTOOL} -g ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -s ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -og ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf 
		#echo "${QCTOOL} -g ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -s ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -og ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf " > ${REGIONALDIR}/${STUDYNAME}_gen2vcf_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
		#qsub -S /bin/bash -N GEN2VCF${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid GENQC${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_gen2vcf_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_gen2vcf_${VARIANT}_excl_${EXCLUSION_TYPE}.log -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_gen2vcf_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
		
		#### Index using Tabix & BGZIP
		#### example: bgzip ${STUDYNAME}_combo_1000g_QC_chr7.7q22.vcf && tabix_v026 -p vcf ${STUDYNAME}_combo_1000g_QC_chr7.7q22.vcf.gz
		echo "* Creating bash-script to submit indexing and gzipping of VCF-file..."
		### for DEBUGGING
		${BGZIP} ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf && ${TABIX} ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz 
		#echo "${BGZIP} ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf && ${TABIX} ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz " > ${REGIONALDIR}/${STUDYNAME}_vcfgz_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
		#qsub -S /bin/bash -N VCFGZ${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid GEN2VCF${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_vcfgz_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_vcfgz_${VARIANT}_excl_${EXCLUSION_TYPE}.log -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_vcfgz_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
	
		echo ""	
		### Running fastQTL
		if [[ ${CHR} -lt 10 ]]; then 
			echo "Processing a variant in region 0${CHR}:${START}-${END}."
			### Running nominal and permutation passes of fastQTL, respectively
			echo "Creating bash-script to submit nominal pass for 'cis-eQTLs'..."
			### for DEBUGGING
			${FASTQTL} --vcf ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${FASTQTLDATA} --region 0${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.log 
			#echo "${FASTQTL} --vcf ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${FASTQTLDATA} --region 0${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.log "> ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			#qsub -S /bin/bash -N QTLnqc${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid VCFGZ${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_NOM} -l h_vmem=${VMEM_NOM} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			echo "Creating bash-script to submit permutation pass for 'cis-eQTLs'..."
			### for DEBUGGING
			${FASTQTL} --vcf ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${FASTQTLDATA} --region 0${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --permute ${PERMSTART} ${PERMEND} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.log
			#echo "${FASTQTL} --vcf ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${FASTQTLDATA} --region 0${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --permute ${PERMSTART} ${PERMEND} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.log "> ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_PERMUTE.sh
			#qsub -S /bin/bash -N QTLpqc${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid VCFGZ${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_PERMUTE.errors -o ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_PERMUTE.output -l h_rt=${QUEUE_PERM} -l h_vmem=${VMEM_PERM} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_PERMUTE.sh
		elif  [[ ${CHR} -ge 10 ]]; then
			echo "Processing a variant in region ${CHR}:${START}-${END}."
			### Running nominal and permutation passes of fastQTL, respectively
			echo "Creating bash-script to submit nominal pass for 'cis-eQTLs'..."
			### for DEBUGGING
			${FASTQTL} --vcf ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${FASTQTLDATA} --region ${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.log 
			#echo "${FASTQTL} --vcf ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${FASTQTLDATA} --region ${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.log "> ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			#qsub -S /bin/bash -N QTLnqc${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid VCFGZ${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_NOM} -l h_vmem=${VMEM_NOM} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			echo ""
			echo "Creating bash-script to submit permutation pass for 'cis-eQTLs'..."
			### for DEBUGGING
			${FASTQTL} --vcf ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${FASTQTLDATA} --region ${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --permute ${PERMSTART} ${PERMEND} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.log 
			#echo "${FASTQTL} --vcf ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${FASTQTLDATA} --region ${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --permute ${PERMSTART} ${PERMEND} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.log "> ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_PERMUTE.sh
			#qsub -S /bin/bash -N QTLpqc${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid VCFGZ${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_PERMUTE.errors -o ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_PERMUTE.output -l h_rt=${QUEUE_PERM} -l h_vmem=${VMEM_PERM} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_PERMUTE.sh
		else
			echo "*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please."	
			exit 1
		fi
		
		echo ""
		
	done < ${REGIONS}

	echo ""
	echo ""
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo "Quality control and parsing of fastQTL results."
	
	while IFS='' read -r REGIONOFINTEREST || [[ -n "$REGIONOFINTEREST" ]]; do
		###	1		2		3	4	5		6		7			8		9
		###Variant	Locus	Chr	BP	BP-1Mb	BP+1Mb	WindowSize	Type	Phenotype
		LINE=${REGIONOFINTEREST}
		VARIANT=$(echo "${LINE}" | awk '{print $1}')
		LOCUS=$(echo "${LINE}" | awk '{print $2}')
		CHR=$(echo "${LINE}" | awk '{print $3}')
		BP=$(echo "${LINE}" | awk '{print $4}')
		START=$(echo "${LINE}" | awk '{print $5}')
		END=$(echo "${LINE}" | awk '{print $6}')
		WINDOWSIZE=$(echo "${LINE}" | awk '{print $7}')
		TYPE=$(echo "${LINE}" | awk '{print $8}')
		PHENOTYPE=$(echo "${LINE}" | awk '{print $9}')
	
		echo "===================================================================="
		echo "Processing ${VARIANT} locus on ${CHR} between ${START} and ${END}..."
		### Make directories for script if they do not exist yet (!!!PREREQUISITE!!!)
		if [ ! -d ${RESULTS}/${VARIANT}_${PROJECTNAME} ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${RESULTS}/${VARIANT}_${PROJECTNAME}
		else
			echo "The regional directory already exists."
		fi
		
		REGIONALDIR=${RESULTS}/${VARIANT}_${PROJECTNAME}
	
		### PERFORMING fastQTL RESULTS QUALITY CONTROL & PARSING
		
		if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
			echo "Creating bash-script to submit 'fastQTL RESULTS QUALITY CONTROL & PARSER v2' on >>> nominal <<< pass results..."
			### FOR DEBUGGING
			Rscript ${FASTQCTLADDON}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t NOM -q MQTL -o ${REGIONALDIR}/ -a ${ORIGINALS}/IlluminaMethylation450K.annotation.txt.gz -j ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats 
			#echo "Rscript ${FASTQCTLADDON}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t NOM -q MQTL -o ${REGIONALDIR}/ -a ${ORIGINALS}/IlluminaMethylation450K.annotation.txt.gz -j ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats  "> ${REGIONALDIR}/${STUDYNAME}_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			#qsub -S /bin/bash -N fastQTLQCnom${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid QTLnqc${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.log -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			echo ""
			echo "Creating bash-script to submit 'fastQTL RESULTS QUALITY CONTROL & PARSER v2' on >>> permutation <<< pass results..."
			### FOR DEBUGGING
			Rscript ${FASTQCTLADDON}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t PERM -q MQTL -o ${REGIONALDIR}/ -a ${ORIGINALS}/IlluminaMethylation450K.annotation.txt.gz -j ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats 
			#echo "Rscript ${FASTQCTLADDON}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t PERM -q MQTL -o ${REGIONALDIR}/ -a ${ORIGINALS}/IlluminaMethylation450K.annotation.txt.gz -j ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats  "> ${REGIONALDIR}/${STUDYNAME}_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			#qsub -S /bin/bash -N fastQTLQCperm${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid QTLpqc${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.log -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			echo ""
		elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
			echo "Creating bash-script to submit 'fastQTL RESULTS QUALITY CONTROL & PARSER v2' on >>> nominal <<< pass results..."
			### FOR DEBUGGING
			Rscript ${FASTQCTLADDON}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t NOM -q EQTL -o ${REGIONALDIR}/ -a ${CTMMEXPRESSIONDATA}/annotation_ctmm_all.csv -j ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats 
			#echo "Rscript ${FASTQCTLADDON}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t NOM -q EQTL -o ${REGIONALDIR}/ -a ${CTMMEXPRESSIONDATA}/annotation_ctmm_all.csv -j ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats  "> ${REGIONALDIR}/${STUDYNAME}_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			#qsub -S /bin/bash -N fastQTLQCnom${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid QTLnqc${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.log -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			echo ""
			echo "Creating bash-script to submit 'fastQTL RESULTS QUALITY CONTROL & PARSER v2' on >>> permutation <<< pass results..."
			### FOR DEBUGGING
			Rscript ${FASTQCTLADDON}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t PERM -q EQTL -o ${REGIONALDIR}/ -a ${CTMMEXPRESSIONDATA}/annotation_ctmm_all.csv -j ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats 
			#echo "Rscript ${FASTQCTLADDON}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t PERM -q EQTL -o ${REGIONALDIR}/ -a ${CTMMEXPRESSIONDATA}/annotation_ctmm_all.csv -j ${REGIONALDIR}/${SNPTESTOUTPUTDATA}_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats  "> ${REGIONALDIR}/${STUDYNAME}_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			#qsub -S /bin/bash -N fastQTLQCperm${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid QTLpqc${STUDYJOBNAME}_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/${STUDYNAME}_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/${STUDYNAME}_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.log -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/${STUDYNAME}_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
			echo ""
		else
			echo "                        *** ERROR *** "
			echo "Something is rotten in the City of Gotham; most likely a typo. "
			echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
			echo "                *** END OF ERROR MESSAGE *** "
			exit 1
		fi
		
	done < ${REGIONS}

	# echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	# echo "Summarising all relevant files into one file."
	# 
	# ### CREATES SUMMARY FILES
	# if [ ! -f ${SUMMARY}/${STUDYNAME}_QC_qtlnom_summary.txt ]; then
	# 	echo "The summary file doesn't exist; Mr. Bourne will make it for you."
	# elif [ ! -f ${SUMMARY}/${STUDYNAME}_QC_qtlperm_summary.txt ]; then
	# 	echo "The summary file doesn't exist; Mr. Bourne will make it for you."
	# else
	# 	echo "The sumary file already exists; Mr. Bourne will re-create it for you."
	# 	rm -v ${SUMMARY}/${STUDYNAME}_QC_qtlnom_summary.txt
	# 	rm -v ${SUMMARY}/${STUDYNAME}_QC_qtlperm_summary.txt
	# fi
	# 
	# if [[ ${STUDY_TYPE} == "AEMS450K1" ]] || [[ ${STUDY_TYPE} == "AEMS450K2" ]]; then
	# 	echo "Making appropriate summary file for results from a mQTL analysis in the '${STUDY_TYPE}'..."
	# 	echo "Locus,ProbeID,VARIANT,Chr,BP,OtherAlleleA,CodedAlleleA,MAF,MAC,CAF,HWE,Info,Imputation,N,Distance_VARIANT_CpG,Chr_CpG,BP_CpG,ProbeType,GeneName,AccessionID_UCSC,GeneGroup_UCSC,CpG_Island_Relation_UCSC,Phantom,DMR,Enhancer,HMM_Island,RegulatoryFeatureName,RegulatoryFeatureGroup,DHS,Beta,SE,Nominal_P,Bonferroni,BenjHoch,Q" > ${SUMMARY}/${STUDYNAME}_QC_qtlnom_summary.txt
	# 	echo "Locus,ProbeID,VARIANT,Chr,BP,OtherAlleleA,CodedAlleleA,MAF,MAC,CAF,HWE,Info,Imputation,N,Distance_VARIANT_CpG,Chr_CpG,BP_CpG,ProbeType,GeneName,AccessionID_UCSC,GeneGroup_UCSC,CpG_Island_Relation_UCSC,Phantom,DMR,Enhancer,HMM_Island,RegulatoryFeatureName,RegulatoryFeatureGroup,DHS,Beta,SE,Nominal_P,Perm_P,ApproxPerm_P,Bonferroni,BenjHoch,Q" > ${SUMMARY}/${STUDYNAME}_QC_qtlperm_summary.txt
	# 	echo ""
	# elif [[ ${STUDY_TYPE} == "CTMM" ]]; then
	# 	echo "Making appropriate summary file for results from an eQTL analysis in the '${STUDY_TYPE}'..."
	# 	echo "Locus,ProbeID,VARIANT,Chr,BP,OtherAlleleA,CodedAlleleA,MAF,MAC,CAF,HWE,Info,Imputation,N,GeneName,EntrezID,Distance_VARIANT_GENE,Chr,GeneTxStart,GeneTxEnd,Beta,SE,Nominal_P,Bonferroni,BenjHoch,Q" > ${SUMMARY}/${STUDYNAME}_QC_qtlnom_summary.txt
	# 	echo "Locus,ProbeID,VARIANT,Chr,BP,OtherAlleleA,CodedAlleleA,MAF,MAC,CAF,HWE,Info,Imputation,N,GeneName,EntrezID,Distance_VARIANT_GENE,Chr,GeneTxStart,GeneTxEnd,Beta,SE,Nominal_P,Perm_P,ApproxPerm_P,Bonferroni,BenjHoch,Q" > ${SUMMARY}/${STUDYNAME}_QC_qtlperm_summary.txt
	# 	echo ""
	# else
	# 	echo "                        *** ERROR *** "
	# 	echo "Something is rotten in the City of Gotham; most likely a typo. "
	# 	echo "Double back, and check you 'study type' please: '${STUDY_TYPE}' does *not* exist."	
	# 	echo "                *** END OF ERROR MESSAGE *** "
	# 	exit 1
	# fi
	# 
	# while IFS='' read -r REGIONOFINTEREST || [[ -n "$REGIONOFINTEREST" ]]; do
	# 	###	1		2		3	4	5		6		7			8		9
	# 	###Variant	Locus	Chr	BP	BP-1Mb	BP+1Mb	WindowSize	Type	Phenotype
	# 	LINE=${REGIONOFINTEREST}
	# 	VARIANT=$(echo "${LINE}" | awk '{print $1}')
	# 	LOCUS=$(echo "${LINE}" | awk '{print $2}')
	# 	CHR=$(echo "${LINE}" | awk '{print $3}')
	# 	BP=$(echo "${LINE}" | awk '{print $4}')
	# 	START=$(echo "${LINE}" | awk '{print $5}')
	# 	END=$(echo "${LINE}" | awk '{print $6}')
	# 	WINDOWSIZE=$(echo "${LINE}" | awk '{print $7}')
	# 	TYPE=$(echo "${LINE}" | awk '{print $8}')
	# 	PHENOTYPE=$(echo "${LINE}" | awk '{print $9}')
	# 
	# 	echo "===================================================================="
	# 	echo "Processing ${VARIANT} locus on ${CHR} between ${START} and ${END}..."
	# 	### Make directories for script if they do not exist yet (!!!PREREQUISITE!!!)
	# 	if [ ! -d ${RESULTS}/${VARIANT}_${PROJECTNAME} ]; then
	# 		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
	# 		mkdir -v ${RESULTS}/${VARIANT}_${PROJECTNAME}
	# 	else
	# 		echo "The regional directory already exists."
	# 	fi
	# 	
	# 	REGIONALDIR=${RESULTS}/${VARIANT}_${PROJECTNAME}
	# 	echo ""
	# 	echo "Copying results to the Summary Directory..."	
	# 	cp -v ${REGIONALDIR}/*_nominal.all.txt ${SUMMARY}/
	# 	cp -v ${REGIONALDIR}/*_perm.P0_05.txt ${SUMMARY}/
	# 	cp -v ${REGIONALDIR}/*.pdf ${SUMMARY}/
	# 	
	# 	echo ""
	# 	echo "Adding all results for the ${VARIANT} locus to the summary file..."
	# 	
	# 	echo ""
	# 	echo "Nominal results..."#tr ',' ' ' | 
	# 	cat ${SUMMARY}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}_nominal.all.txt | tail -n +2 | awk -v LOCUS_VARIANT=$VARIANT '{ print LOCUS_VARIANT, $0 }' OFS=","  >> ${SUMMARY}/${STUDYNAME}_QC_qtlnom_summary.txt
	# 	gzip -v ${SUMMARY}/${STUDYNAME}_QC_qtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}_nominal.all.txt
	# 	echo ""
	# 	echo "Permutation results..."
	# 	cat ${SUMMARY}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}_perm.P0_05.txt | tail -n +2 | awk -v LOCUS_VARIANT=$VARIANT '{ print LOCUS_VARIANT, $0 }' OFS=","  >> ${SUMMARY}/${STUDYNAME}_QC_qtlperm_summary.txt
	# 	gzip -v ${SUMMARY}/${STUDYNAME}_QC_qtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}_perm.P0_05.txt
	# 	
	# done < ${REGIONS}

	###ZIPPING FINAL SUMMARY RESULTS
	# gzip -v ${SUMMARY}/${STUDYNAME}_QC_qtlnom_summary.txt
	# gzip -v ${SUMMARY}/${STUDYNAME}_QC_qtlperm_summary.txt

	# echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
	# echo "Parse nominal results to get Hits per Locus and (mapped) Gene, and input files"
	# echo "for LocusZoom v1.2+."
	# 
	# ### First we will collect all the nominal association results.
	# echo ""
	# echo "Parsing nominal results..."
	# cd ${SUMMARY}
	# pwd
	# ###module load python
	# ###python ${FASTQTLPARSER} ${SUMMARY}/${STUDYNAME}_QC_qtlnom_summary.txt.gz

	### Now we will start plotting per locus each gene-probe-pair.
	#echo ""
	#
	# while IFS='' read -r REGIONOFINTEREST || [[ -n "$REGIONOFINTEREST" ]]; do
	# ### 1		2		3	4	5		6		7			8		9
	# ### Variant	Locus	Chr	BP	BP-1Mb	BP+1Mb	WindowSize	Type	Phenotype
	# LINE=${REGIONOFINTEREST}
	# VARIANT=$(echo "${LINE}" | awk '{print $1}')
	# LOCUS=$(echo "${LINE}" | awk '{print $2}')
	# CHR=$(echo "${LINE}" | awk '{print $3}')
	# BP=$(echo "${LINE}" | awk '{print $4}')
	# START=$(echo "${LINE}" | awk '{print $5}')
	# END=$(echo "${LINE}" | awk '{print $6}')
	# WINDOWSIZE=$(echo "${LINE}" | awk '{print $7}')
	# TYPE=$(echo "${LINE}" | awk '{print $8}')
	# PHENOTYPE=$(echo "${LINE}" | awk '{print $9}')
	# 
	# echo "===================================================================="
	# echo "        INITIALISING LOCUSZOOM PLOTTING FOR ${VARIANT} LOCUS"
	# echo "===================================================================="
	# echo ""
	# 
	# ### Getting only the top part of the variant-gene-probeid list
	# cat ${SUMMARY}/_loci/${VARIANT}.txt | tail -n +2 > ${SUMMARY}/_loci/${VARIANT}.LZ.txt
	# LOCUSHITS=${SUMMARY}/_loci/${VARIANT}.LZ.txt
	# echo "These are the hits we're interested in for ${VARIANT}..."
	# echo ""
	# cat ${LOCUSHITS}
	# 
	# echo ""
	# while IFS='' read -r VARIANTGENEPROBE || [[ -n "$VARIANTGENEPROBE" ]]; do
	# 	### 1		2			3		4			5
	# 	### Locus	GeneName	ProbeID	N_Variants	N_Significant
	# 	LINE=${VARIANTGENEPROBE}
	# 	LOCUSVARIANT=$(echo "${LINE}" | awk '{print $1}')
	# 	GENENAME=$(echo "${LINE}" | awk '{print $2}')
	# 	PROBEID=$(echo "${LINE}" | awk '{print $3}')
	# 	N_VARIANTS=$(echo "${LINE}" | awk '{print $4}')
	# 	N_SIGNIFICANT=$(echo "${LINE}" | awk '{print $5}')
	# 
	# 	echo "===================================================================="
	# 	echo "Plotting results for the ${LOCUSVARIANT} locus on ${CHR}:${START}-${END}."
	# 	echo "	* Plotting association results for ${GENENAME} and ${PROBEID}."
	# 	echo "	* Total number of variants analysed: 	${N_VARIANTS}."
	# 	echo "	* Total number of significant variants:	${N_SIGNIFICANT}."
	# 	
	# 	echo ""
	# 	### FOR DEBUGGING
	# 	### echo "Checking existence of proper input file..."
	# 	### ls -lh ${SUMMARY}/_probes/${LOCUSVARIANT}_${GENENAME}_${PROBEID}.lz
	# 	### echo ""
	# 	### echo "Head:"
	# 	### head ${SUMMARY}/_probes/${LOCUSVARIANT}_${GENENAME}_${PROBEID}.lz
	# 	### echo ""
	# 	### echo "Tail:"
	# 	### tail ${SUMMARY}/_probes/${LOCUSVARIANT}_${GENENAME}_${PROBEID}.lz
	# 	### echo ""
	# 	### echo "Row count:"
	# 	### cat ${SUMMARY}/_probes/${LOCUSVARIANT}_${GENENAME}_${PROBEID}.lz | wc -l
	# 	
	# 	### Setting up LocusZoom v1.2+ plotting
	# 	### Some general settings
	# 	LOCUSZOOM_SETTINGS="ldColors=\"#595A5C,#4C81BF,#1396D8,#C5D220,#F59D10,red,#9A3480\" showRecomb=TRUE dCol='r^2' drawMarkerNames=FALSE refsnpTextSize=0.8 geneFontSize=0.6 showRug=FALSE showAnnot=FALSE showRefsnpAnnot=TRUE showGenes=TRUE clean=TRUE bigDiamond=TRUE ymax=12 rfrows=10 refsnpLineWidth=2"
	# 	### The proper genome-build
	# 	LDMAP="--pop EUR --build hg19 --source 1000G_March2012"
	# 	### Directory prefix
	# 	PREFIX="${LOCUSVARIANT}_${GENE}_${PROBEID}_excl_${EXCLUSION_TYPE}_"
	# 	
	# 	### Actual plotting
	# 	${LZ13} --metal ${SUMMARY}/_probes/${LOCUSVARIANT}_${GENE}_${PROBEID}.lz --markercol MarkerName --pvalcol P-value --delim tab --chr ${CHR} --start ${START} --end ${END} ${LDMAP} ${LOCUSZOOM_SETTINGS} --prefix=${PREFIX} theme=publication title="${LOCUSVARIANT} - ${GENENAME} (${PROBEID})" 
	# 
	# done < ${LOCUSHITS}
	# 
	# ### rm -v ${SUMMARY}/_loci/${VARIANT}.LZ.txt
	# 
	# done < ${REGIONS}


	### CODE BELOW PROBABLY OBSOLETE
	###done ### COMMENT OUT IF YOU DON'T WANT TO ANALYZE ALL STRATIFICATION OPTIONS
	### CODE ABOVE PROBABLY OBSOLETE

### END of if-else statement for the number of command-line arguments passed ###
fi

script_copyright_message
	

	### EXCLUSION LISTS NOTES
	### Please note that it is okay to have more (or different) samples in the *genotype* data 
	### as compared to the *phenotype* data. However, it is *NOT* okay to have more 
	### (or different) samples in the *phenotype* data as compared to the *genotype* data!!!
	### In other words: remove from the BED files - *BEFORE* while making them!!! - the
	### samples that do *NOT* have genotype data!!!
	### FOR DEBUGGING USEFULL:
	###EXCLUSION_NONAEMS450K1="${GENETICDATA}/exclude_nonAEMS450K1.list"
	###EXCLUSION_NONAEMS450K2="${GENETICDATA}/exclude_nonAEMS450K2.list"


	### EXCLUSION LISTS for fastQTL
	###
	### CTMM
	### - exclusion_nonCTMM_forFastQTL.list
	### - exclusion_nonCTMM_FEMALES_forFastQTL.list
	### - exclusion_nonCTMM_MALES_forFastQTL.list
	### - exclusion_nonCTMM_SMOKER_forFastQTL.list
	### - exclusion_nonCTMM_nonSMOKER_forFastQTL.list
	### - exclusion_nonCTMM_T2D_forFastQTL.list
	### - exclusion_nonCTMM_nonT2D_forFastQTL.list
	###
	### AEGS
	### For AEGS the exclusion lists have exactly the same identifiers.

	### EXCLUSION LISTS for SNPTEST
	###
	### CTMM
	### - exclusion_nonCTMM.list
	### - exclusion_nonCTMM_FEMALES.list
	### - exclusion_nonCTMM_MALES.list
	### - exclusion_nonCTMM_SMOKER.list
	### - exclusion_nonCTMM_nonSMOKER.list
	### - exclusion_nonCTMM_T2D.list
	### - exclusion_nonCTMM_nonT2D.list
	###
	### AEGS
	### - exclusion_Females.list
	### - exclusion_Males.list
	### - exclusion_nonAEGS.list
	### - exclusion_nonCEA_AEGS.list
	### - exclusion_nonCEA.list
	### - exclusion_nonCEA_Females.list
	### - exclusion_nonCEA_Males.list
	### - exclusion_nonCEA_SMOKER.list
	### - exclusion_nonCEA_nonSMOKER.list
	### - exclusion_nonCEA_T2D.list
	### - exclusion_nonCEA_nonT2D.list
	### -- AEGS SPECIFIC --
	### - exclusion_nonCEA_CKD.list
	### - exclusion_nonCEA_nonCKD.list
	### - exclusion_nonCEA_post2007.list
	### - exclusion_nonCEA_pre2007.list
