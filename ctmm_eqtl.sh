#!/bin/bash
#
# You can use the variables below (indicated by "#$") to set some things for the 
# submission system.
#$ -S /bin/bash # the type of BASH you'd like to use
#$ -N ctmm_eqtl_v1_20151223 # the name of this script
#$ -hold_jid some_other_basic_bash_script # the current script (basic_bash_script) will hold until some_other_basic_bash_script has finished
#$ -o /hpc/dhl_ec/svanderlaan/projects/ctmm/ctmm_eqtl.log # the log file of this job
#$ -e /hpc/dhl_ec/svanderlaan/projects/ctmm/ctmm_eqtl.errors # the error file of this job
#$ -l h_rt=24:00:00 # h_rt=[max time, hh:mm:ss, e.g. 02:02:01] - this is the time you think the script will take
#$ -l h_vmem=16G #  h_vmem=[max. mem, e.g. 45G] - this is the amount of memory you think your script will use
#$ -l tmpspace=32G # this is the amount of temporary space you think your script will use
#$ -M s.w.vanderlaan-2@umcutrecht.nl # you can send yourself emails when the job is done; "-M" and "-m" go hand in hand
#$ -m beas # you can choose: Emails can be sent upon start ( b) and/or end ( e) and/or abortion ( a) and/or suspension ( s) of your job. 
#$ -cwd # set the job start to the current directory - so all the things in this script are relative to the current directory!!!
#
# Another useful tip: you can set a job to run after another has finished. Name the job 
# with "-N SOMENAME" and hold the other job with -hold_jid SOMENAME". 
# Further instructions: https://wiki.bioinformatics.umcutrecht.nl/bin/view/HPC/HowToS#Run_a_job_after_your_other_jobs
#
# The command 'clear' cleares the screen.

# It is good practice to properly name and annotate your script for future reference for
# yourself and others. Trust me, you'll forget why and how you made this!!!
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "                  EXPRESSION QUANTITATIVE TRAIT LOCUS ANALYSIS OF CTMM"
echo "                                 version 1.0 (20160308)"
echo ""
echo "* Written by  : Sander W. van der Laan"
echo "* E-mail      : s.w.vanderlaan-2@umcutrecht.nl"
echo "* Last update : 2016-03-08"
echo "* Version     : ctmm_eqtl_v1_20160308"
echo ""
echo "* Description : This script will set some directories, execute something in a for "
echo "                loop, and will then submit this in a job."
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Today's: "$(date)
TODAY=$(date +"%Y%m%d")
echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "The following directories are set."
SOFTWARE=/hpc/local/CentOS7/dhl_ec/software
QCTOOL=${SOFTWARE}/qctool_v1.5-linux-x86_64-static/qctool
SNPTEST252=${SOFTWARE}/snptest_v2.5.2_CentOS6.5_x86_64_static/snptest_v2.5.2
FASTQTL=${SOFTWARE}/fastqtl_v2.184
FASTQTLPARSER=${SOFTWARE}/fastqtl_supportives/NominalResultsParser.py
LZ13=${SOFTWARE}/locuszoom_1.3/bin/locuszoom
BGZIP=${SOFTWARE}/htslib-1.3/bgzip
TABIX=${SOFTWARE}/htslib-1.3/tabix

### CTMM SPECIFIC
ORIGINALS=/hpc/dhl_ec/data/_ctmm_originals
GENETICDATA=${ORIGINALS}/CTMMAxiomTX_IMPUTE2_1000Gp3_GoNL5
EXPRESSIONDATA=${ORIGINALS}/CTMMHumanHT12v4r2_15002873B
PROJECT=/hpc/dhl_ec/svanderlaan/projects/ctmm

echo ""
echo "Software directory                                                   ${SOFTWARE}"
echo "Where \"qctool\" resides                                               ${QCTOOL}"
echo "Where \"fastQTL\" resides                                              ${FASTQTL}"
echo "Where \"bgzip\" resides                                                ${BGZIP}"
echo "Where \"tabix\" resides                                                ${TABIX}"
echo "Where \"snptest 2.5.2\" resides                                        ${SNPTEST252}"
echo ""
echo "Original CTMM data directory                                         ${ORIGINALS}"
echo "CTMM genetic data directory                                          ${GENETICDATA}"
echo "CTMM expression data directory                                       ${EXPRESSIONDATA}"
echo "Project directory                                                    ${PROJECT}"
echo ""
echo "We will run this script on                                           ${TODAY}"
echo ""
echo "Additional fastQTL specific settings."

### FASTQTL & QCTOOL SETTINGS
SEEDNO=91149216
REGIONS=${PROJECT}/regions_for_eqtl.txt # regions_for_eqtl.txt
PERMSTART="1000"
PERMEND="1000000"
MAF="0.005"
INFO="0.9"
HWE="6"

### COVARIATES FILE
COVARIATES="${GENETICDATA}/covariates_ctmm.all.cov"

### EXCLUSION LISTS for fastQTL
### - exclusion_nonCTMM_forFastQTL.list
### - exclusion_nonCTMM_FEMALES_forFastQTL.list
### - exclusion_nonCTMM_MALES_forFastQTL.list
### - exclusion_nonCTMM_SMOKER_forFastQTL.list
### - exclusion_nonCTMM_nonSMOKER_forFastQTL.list
### - exclusion_nonCTMM_T2D_forFastQTL.list
### - exclusion_nonCTMM_nonT2D_forFastQTL.list

### EXCLUSION LISTS for SNPTEST
### - exclusion_nonCTMM.list
### - exclusion_nonCTMM_FEMALES.list
### - exclusion_nonCTMM_MALES.list
### - exclusion_nonCTMM_SMOKER.list
### - exclusion_nonCTMM_nonSMOKER.list
### - exclusion_nonCTMM_T2D.list
### - exclusion_nonCTMM_nonT2D.list

#for TYPE in DEFAULT MALES FEMALES SMOKER NONSMOKER T2D NONT2D; do

EXCLUSION_TYPE="DEFAULT" # DEFAULT/SMOKER/NONSMOKER/MALES/FEMALES/T2D/NONT2D

### QSUB SETTINGS
QUEUE_QCTOOL="01:00:00"
VMEM_QCTOOL="4G"
QUEUE_NOM="08:00:00"
VMEM_NOM="4G"
QUEUE_PERM="24:00:00"
VMEM_PERM="4G"

EMAIL="s.w.vanderlaan-2@umcutrecht.nl"
MAILTYPE="as"

PROJECTNAME="CAD"

echo ""
echo "Seed number                                                           ${SEEDNO}"

### AEMS - CTMM SPECIFIC
### CREATES RESULTS AND SUMMARY DIRECTORIES; SETS EXCLUSION-FILES
### DEFAULT
if [[ ${EXCLUSION_TYPE} == "DEFAULT" ]]; then
	# Results directory
	if [ ! -d ${PROJECT}/eqtl ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/eqtl
		else
			echo "The regional directory already exists."
		fi
	RESULTS=${PROJECT}/eqtl
	# Summary directory
	if [ ! -d ${PROJECT}/eqtl_summary ]; then
		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
		mkdir -v ${PROJECT}/eqtl_summary
	else
		echo "The regional directory already exists."
	fi
	SUMMARY=${PROJECT}/eqtl_summary
	# Exclusion files
	EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM.list"
	EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_forFastQTL.list"
	EXCLUSION_COV="${PROJECT}/excl_cov.txt"

### MALES ONLY ANALYSIS
elif [[ ${EXCLUSION_TYPE} == "FEMALES" ]]; then
	# Results directory
	if [ ! -d ${PROJECT}/eqtl_males ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/eqtl_males
		else
			echo "The regional directory already exists."
		fi
	RESULTS=${PROJECT}/eqtl_males
	# Summary directory
	if [ ! -d ${PROJECT}/eqtl_summary_males ]; then
		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
		mkdir -v ${PROJECT}/eqtl_summary_males
	else
		echo "The regional directory already exists."
	fi
	SUMMARY=${PROJECT}/eqtl_summary_males
	# Exclusion files
	EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_FEMALES.list"
	EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_FEMALES_forFastQTL.list"
	EXCLUSION_COV="${PROJECT}/excl_cov_sex.txt"
	
### FEMALES ONLY ANALYSIS
elif [[ ${EXCLUSION_TYPE} == "MALES" ]]; then		
	# Results directory
	if [ ! -d ${PROJECT}/eqtl_females ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/eqtl_females
		else
			echo "The regional directory already exists."
		fi
	RESULTS=${PROJECT}/eqtl_females
	# Summary directory	
	if [ ! -d ${PROJECT}/eqtl_summary_females ]; then
		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
		mkdir -v ${PROJECT}/eqtl_summary_females
	else
		echo "The regional directory already exists."
	fi
	SUMMARY=${PROJECT}/eqtl_summary_females
	# Exclusion files
	EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_MALES.list"
	EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_MALES_forFastQTL.list"
	EXCLUSION_COV="${PROJECT}/excl_cov_sex.txt"

### SMOKER ONLY ANALYSIS
elif [[ ${EXCLUSION_TYPE} == "NONSMOKER" ]]; then		
	# Results directory
	if [ ! -d ${PROJECT}/eqtl_smoker ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/eqtl_smoker
		else
			echo "The regional directory already exists."
		fi
	RESULTS=${PROJECT}/eqtl_smoker
	# Summary directory	
	if [ ! -d ${PROJECT}/eqtl_summary_smoker ]; then
		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
		mkdir -v ${PROJECT}/eqtl_summary_smoker
	else
		echo "The regional directory already exists."
	fi
	SUMMARY=${PROJECT}/eqtl_summary_smoker
	# Exclusion files
	EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_nonSMOKER.list"
	EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_nonSMOKER_forFastQTL.list"
	EXCLUSION_COV="${PROJECT}/excl_cov.txt"
	
### NONSMOKER ONLY ANALYSIS
elif [[ ${EXCLUSION_TYPE} == "SMOKER" ]]; then		
	# Results directory
	if [ ! -d ${PROJECT}/eqtl_nonsmoker ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/eqtl_nonsmoker
		else
			echo "The regional directory already exists."
		fi
	RESULTS=${PROJECT}/eqtl_nonsmoker
	# Summary directory	
	if [ ! -d ${PROJECT}/eqtl_summary_nonsmoker ]; then
		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
		mkdir -v ${PROJECT}/eqtl_summary_nonsmoker
	else
		echo "The regional directory already exists."
	fi
	SUMMARY=${PROJECT}/eqtl_summary_nonsmoker
	# Exclusion files
	EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_SMOKER.list"
	EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_SMOKER_forFastQTL.list"
	EXCLUSION_COV="${PROJECT}/excl_cov.txt"
	
### NON-TYPE 2 DIABETES ONLY ANALYSIS
elif [[ ${EXCLUSION_TYPE} == "T2D" ]]; then		
	# Results directory
	if [ ! -d ${PROJECT}/eqtl_nont2d ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/eqtl_nont2d
		else
			echo "The regional directory already exists."
		fi
	RESULTS=${PROJECT}/eqtl_nont2d
	# Summary directory	
	if [ ! -d ${PROJECT}/eqtl_summary_nont2d ]; then
		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
		mkdir -v ${PROJECT}/eqtl_summary_nont2d
	else
		echo "The regional directory already exists."
	fi
	SUMMARY=${PROJECT}/eqtl_summary_nont2d
	# Exclusion files
	EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_T2D.list"
	EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_T2D_forFastQTL.list"
	EXCLUSION_COV="${PROJECT}/excl_cov.txt"
	
### TYPE 2 DIABETES ONLY ANALYSIS
elif [[ ${EXCLUSION_TYPE} == "NONT2D" ]]; then		
	# Results directory
	if [ ! -d ${PROJECT}/eqtl_t2d ]; then
			echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
			mkdir -v ${PROJECT}/eqtl_t2d
		else
			echo "The regional directory already exists."
		fi
	RESULTS=${PROJECT}/eqtl_t2d
	# Summary directory	
	if [ ! -d ${PROJECT}/eqtl_summary_t2d ]; then
		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
		mkdir -v ${PROJECT}/eqtl_summary_t2d
	else
		echo "The regional directory already exists."
	fi
	SUMMARY=${PROJECT}/eqtl_summary_t2d
	# Exclusion files
	EXCLUSION_SNPTEST="${GENETICDATA}/exclusion_nonCTMM_nonT2D.list"
	EXCLUSION_FASTQTL="${GENETICDATA}/exclusion_nonCTMM_nonT2D_forFastQTL.list"
	EXCLUSION_COV="${PROJECT}/excl_cov.txt"

else
	echo "Something's rotten in the City of Gotham. Please, double back."
fi


echo ""
echo "The list of regions to investigate                                    "
	echo "* Variant	Locus	Chr	BP	StartRange	EndRange	WindowSize	Type	Phenotype"
	while IFS='' read -r REGION || [[ -n "$REGION" ]]; do
	LINE=${REGION}
	echo "* ${LINE}"
	done < ${REGIONS}
echo ""

#echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
#echo "Extracting CARDIoGRAMplusC4D (2015) loci with a 1,000,000 bp range."
#while IFS='' read -r REGIONOFINTEREST || [[ -n "$REGIONOFINTEREST" ]]; do
#	#	1		2		3	4	5		6		7			8		9
#	#Variant	Locus	Chr	BP	StartRange	EndRange	WindowSize	Type	Phenotype
#	LINE=${REGIONOFINTEREST}
#	VARIANT=$(echo "${LINE}" | awk '{print $1}')
#	LOCUS=$(echo "${LINE}" | awk '{print $2}')
#	CHR=$(echo "${LINE}" | awk '{print $3}')
#	BP=$(echo "${LINE}" | awk '{print $4}')
#	START=$(echo "${LINE}" | awk '{print $5}')
#	END=$(echo "${LINE}" | awk '{print $6}')
#	WINDOWSIZE=$(echo "${LINE}" | awk '{print $7}')
#	TYPE=$(echo "${LINE}" | awk '{print $8}')
#	PHENOTYPE=$(echo "${LINE}" | awk '{print $9}')
#
#	echo "===================================================================="
#	echo "Processing ${VARIANT} locus on ${CHR} between ${START} and ${END}..."
#	### Make directories for script if they do not exist yet (!!!PREREQUISITE!!!)
#	if [ ! -d ${RESULTS}/${VARIANT}_${PROJECTNAME} ]; then
#		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
#		mkdir -v ${RESULTS}/${VARIANT}_${PROJECTNAME}
#	else
#		echo "The regional directory already exists."
#	fi
#	
#	REGIONALDIR=${RESULTS}/${VARIANT}_${PROJECTNAME}
#	
#	### Extraction relevant regions for CTMM eQTL analysis using fastQTL
#	echo "Creating bash-script to submit qctool extraction of region ${CHR}:${START}-${END} near ${LOCUS}..."
#	echo "${QCTOOL} -g ${GENETICDATA}/ctmm_1kGp3GoNL5_RAW_chr${CHR}.bgen -s ${GENETICDATA}/ctmm_1kGp3GoNL5_RAW_chr${CHR}.sample -og ${REGIONALDIR}/ctmm_1kGp3GoNL5_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -os ${REGIONALDIR}/ctmm_1kGp3GoNL5_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -incl-range ${CHR}:${START}-${END} "> ${REGIONALDIR}/ctmm_genex_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	qsub -S /bin/bash -N GENEX_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_genex_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/ctmm_genex_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_genex_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	
#	echo "Creating bash-script to submit qctool filtering of the ${LOCUS} based on MAF > ${MAF}, INFO > ${INFO} and HWE -log10(p) > ${HWE}..."
#	echo "${QCTOOL} -g ${REGIONALDIR}/ctmm_1kGp3GoNL5_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -og ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -maf ${MAF} 1 -info ${INFO} 1 -hwe ${HWE} "> ${REGIONALDIR}/ctmm_genqc_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	qsub -S /bin/bash -N GENQC_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid GENEX_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_genqc_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/ctmm_genqc_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_genqc_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	
#	### Calculating statistics CTMM
#	echo "Creating bash-script to submit to calculate summary statistics of region ${CHR}:${START}-${END}..."
#	echo "${SNPTEST252} -data ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz ${GENETICDATA}/ctmm_phenocov.sample -summary_stats_only -hwe -exclude_samples ${EXCLUSION_SNPTEST} -o ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats "> ${REGIONALDIR}/ctmm_genstats_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	qsub -S /bin/bash -N GENSTAT_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid GENQC_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_genstats_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/ctmm_genstats_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_genstats_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#
#	### Make VCF of CTMM
#	### example: qctool_v15 -g ctmm_1kGp3GoNL5_QC_chr7.7q22.gen.gz -s ctmm_phenocov.sample -og ctmm_1kGp3GoNL5_QC_chr7.7q22.vcf
#	echo "Creating bash-script to submit VCF-file generation..."
#	echo "${QCTOOL} -g ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.gen.gz -s ${REGIONALDIR}/ctmm_1kGp3GoNL5_RAW_${VARIANT}_excl_${EXCLUSION_TYPE}.sample -og ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf " > ${REGIONALDIR}/ctmm_gen2vcf_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	qsub -S /bin/bash -N GEN2VCF_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid GENQC_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_gen2vcf_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/ctmm_gen2vcf_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_gen2vcf_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	
#	### Index using Tabix & BGZIP
#	### example: bgzip aegs_combo_1000g_QC_chr7.7q22.vcf && tabix_v026 -p vcf aegs_combo_1000g_QC_chr7.7q22.vcf.gz
#	echo "Creating bash-script to submit indexing and gzipping of VCF-file..."
#	echo "${BGZIP} ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf && ${TABIX} ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz " > ${REGIONALDIR}/ctmm_vcfgz_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	qsub -S /bin/bash -N VCFGZ_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid GEN2VCF_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_vcfgz_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/ctmm_vcfgz_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_vcfgz_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#		
#	### Running fastQTL
#	if [[ ${CHR} -lt 10 ]]; then 
#		echo "Processing a variant in region 0${CHR}:${START}-${END}."
#		### Running nominal and permutation passes of fastQTL, respectively
#		echo "Creating bash-script to submit nominal pass for 'cis-eQTLs'..."
#		echo "${FASTQTL} --vcf ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${EXPRESSIONDATA}/phenotype_ctmm.all.bed.gz --region 0${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.log "> ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#		qsub -S /bin/bash -N eQTLnqc_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid VCFGZ_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_NOM} -l h_vmem=${VMEM_NOM} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#		echo ""
#		echo "Creating bash-script to submit permutation pass for 'cis-eQTLs'..."
#		echo "${FASTQTL} --vcf ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${EXPRESSIONDATA}/phenotype_ctmm.all.bed.gz --region 0${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --permute ${PERMSTART} ${PERMEND} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.log "> ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_PERMUTE.sh
#		qsub -S /bin/bash -N eQTLpqc_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid VCFGZ_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_PERMUTE.errors -o ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_PERMUTE.output -l h_rt=${QUEUE_PERM} -l h_vmem=${VMEM_PERM} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_PERMUTE.sh
#	elif  [[ ${CHR} -ge 10 ]]; then
#		echo "Processing a variant in region ${CHR}:${START}-${END}."
#		### Running nominal and permutation passes of fastQTL, respectively
#		echo "Creating bash-script to submit nominal pass for 'cis-eQTLs'..."
#		echo "${FASTQTL} --vcf ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${EXPRESSIONDATA}/phenotype_ctmm.all.bed.gz --region ${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.log "> ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#		qsub -S /bin/bash -N eQTLnqc_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid VCFGZ_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_NOM} -l h_vmem=${VMEM_NOM} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#		echo ""
#		echo "Creating bash-script to submit permutation pass for 'cis-eQTLs'..."
#		echo "${FASTQTL} --vcf ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.vcf.gz --bed ${EXPRESSIONDATA}/phenotype_ctmm.all.bed.gz --region ${CHR}:${START}-${END} --seed ${SEEDNO} --window ${WINDOWSIZE} --permute ${PERMSTART} ${PERMEND} --exclude-samples ${EXCLUSION_FASTQTL} --exclude-covariates ${EXCLUSION_COV} --cov ${COVARIATES} --out ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz --log ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.log "> ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_PERMUTE.sh
#		qsub -S /bin/bash -N eQTLpqc_${VARIANT}_excl_${EXCLUSION_TYPE} -hold_jid VCFGZ_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_PERMUTE.errors -o ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_PERMUTE.output -l h_rt=${QUEUE_PERM} -l h_vmem=${VMEM_PERM} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_PERMUTE.sh
#	else
#		echo "*** ERROR *** Something is rotten in the City of Gotham; most likely a typo. Double back, please."	
#		exit 1
#	fi
#	echo ""
#	
#done < ${REGIONS}

#echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
#echo "Quality control and parsing of fastQTL results."
#
#while IFS='' read -r REGIONOFINTEREST || [[ -n "$REGIONOFINTEREST" ]]; do
#	#	1		2		3	4	5		6		7			8		9
#	#Variant	Locus	Chr	BP	BP-1Mb	BP+1Mb	WindowSize	Type	Phenotype
#	LINE=${REGIONOFINTEREST}
#	VARIANT=$(echo "${LINE}" | awk '{print $1}')
#	LOCUS=$(echo "${LINE}" | awk '{print $2}')
#	CHR=$(echo "${LINE}" | awk '{print $3}')
#	BP=$(echo "${LINE}" | awk '{print $4}')
#	START=$(echo "${LINE}" | awk '{print $5}')
#	END=$(echo "${LINE}" | awk '{print $6}')
#	WINDOWSIZE=$(echo "${LINE}" | awk '{print $7}')
#	TYPE=$(echo "${LINE}" | awk '{print $8}')
#	PHENOTYPE=$(echo "${LINE}" | awk '{print $9}')
#
#	echo "===================================================================="
#	echo "Processing ${VARIANT} locus on ${CHR} between ${START} and ${END}..."
#	### Make directories for script if they do not exist yet (!!!PREREQUISITE!!!)
#	if [ ! -d ${RESULTS}/${VARIANT}_${PROJECTNAME} ]; then
#		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
#		mkdir -v ${RESULTS}/${VARIANT}_${PROJECTNAME}
#	else
#		echo "The regional directory already exists."
#	fi
#	
#	REGIONALDIR=${RESULTS}/${VARIANT}_${PROJECTNAME}
#
#	### PERFORMING fastQTL RESULTS QUALITY CONTROL & PARSING
#	echo "Creating bash-script to submit 'fastQTL RESULTS QUALITY CONTROL & PARSER v1' on nominal pass results..."
#	#Rscript ${PROJECT}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t NOM -q EQTL -o ${REGIONALDIR}/ -a ${EXPRESSIONDATA}/annotation_ctmm_all.csv -j ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats 
#	echo "Rscript ${PROJECT}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t NOM -q EQTL -o ${REGIONALDIR}/ -a ${EXPRESSIONDATA}/annotation_ctmm_all.csv -j ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats  "> ${REGIONALDIR}/ctmm_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	qsub -S /bin/bash -N fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/ctmm_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_fastQTLQCnom_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	echo ""
#	echo "Creating bash-script to submit 'fastQTL RESULTS QUALITY CONTROL & PARSER v1' on permutation pass results..."
#	#Rscript ${PROJECT}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t PERM -q EQTL -o ${REGIONALDIR}/ -a ${EXPRESSIONDATA}/annotation_ctmm_all.csv -j ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats 
#	echo "Rscript ${PROJECT}/fastQTL_QC.R -p ${PROJECT} -r ${REGIONALDIR}/ctmm_QC_eqtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}.txt.gz -t PERM -q EQTL -o ${REGIONALDIR}/ -a ${EXPRESSIONDATA}/annotation_ctmm_all.csv -j ${REGIONALDIR}/ctmm_1kGp3GoNL5_QC_${VARIANT}_excl_${EXCLUSION_TYPE}.stats  "> ${REGIONALDIR}/ctmm_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#	qsub -S /bin/bash -N fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE} -e ${REGIONALDIR}/ctmm_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.errors -o ${REGIONALDIR}/ctmm_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.output -l h_rt=${QUEUE_QCTOOL} -l h_vmem=${VMEM_QCTOOL} -M ${EMAIL} -m ${MAILTYPE} -wd ${REGIONALDIR} ${REGIONALDIR}/ctmm_fastQTLQCperm_${VARIANT}_excl_${EXCLUSION_TYPE}.sh
#
#done < ${REGIONS}

#echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
#echo "Summarising all relevant files into one file."
#
#### CREATES SUMMARY FILES
#if [ ! -f ${SUMMARY}/ctmm_QC_eqtlnom_summary.txt ]; then
#	echo "The summary file doesn't exist; Mr. Bourne will make it for you."
#elif [ ! -f ${SUMMARY}/ctmm_QC_eqtlperm_summary.txt ]; then
#	echo "The summary file doesn't exist; Mr. Bourne will make it for you."
#else
#	echo "The sumary file already exists; Mr. Bourne will re-create it for you."
#	rm -v ${SUMMARY}/ctmm_QC_eqtlnom_summary.txt
#	rm -v ${SUMMARY}/ctmm_QC_eqtlperm_summary.txt
#fi
#echo "Locus	ProbeID	VARIANT	Chr	BP	OtherAlleleA	CodedAlleleA	MAF	MAC	CAF	HWE	Info	Imputation	N	GeneName	EntrezID	Distance_VARIANT_GENE	Chr	GeneTxStart	GeneTxEnd	Beta	SE	Nominal_P	Bonferroni	BenjHoch	Q" > ${SUMMARY}/ctmm_QC_eqtlnom_summary.txt
#echo "Locus	ProbeID	VARIANT	Chr	BP	OtherAlleleA	CodedAlleleA	MAF	MAC	CAF	HWE	Info	Imputation	N	GeneName	EntrezID	Distance_VARIANT_GENE	Chr	GeneTxStart	GeneTxEnd	Beta	SE	Nominal_P	Perm_P	ApproxPerm_P	Bonferroni	BenjHoch	Q" > ${SUMMARY}/ctmm_QC_eqtlperm_summary.txt
#
#while IFS='' read -r REGIONOFINTEREST || [[ -n "$REGIONOFINTEREST" ]]; do
#	#	1		2		3	4	5		6		7			8		9
#	#Variant	Locus	Chr	BP	BP-1Mb	BP+1Mb	WindowSize	Type	Phenotype
#	LINE=${REGIONOFINTEREST}
#	VARIANT=$(echo "${LINE}" | awk '{print $1}')
#	LOCUS=$(echo "${LINE}" | awk '{print $2}')
#	CHR=$(echo "${LINE}" | awk '{print $3}')
#	BP=$(echo "${LINE}" | awk '{print $4}')
#	START=$(echo "${LINE}" | awk '{print $5}')
#	END=$(echo "${LINE}" | awk '{print $6}')
#	WINDOWSIZE=$(echo "${LINE}" | awk '{print $7}')
#	TYPE=$(echo "${LINE}" | awk '{print $8}')
#	PHENOTYPE=$(echo "${LINE}" | awk '{print $9}')
#
#	echo "===================================================================="
#	echo "Processing ${VARIANT} locus on ${CHR} between ${START} and ${END}..."
#	### Make directories for script if they do not exist yet (!!!PREREQUISITE!!!)
#	if [ ! -d ${RESULTS}/${VARIANT}_${PROJECTNAME} ]; then
#		echo "The regional directory doesn't exist; Mr. Bourne will make it for you."
#		mkdir -v ${RESULTS}/${VARIANT}_${PROJECTNAME}
#	else
#		echo "The regional directory already exists."
#	fi
#	
#	REGIONALDIR=${RESULTS}/${VARIANT}_${PROJECTNAME}
#	echo ""
#	echo "Copying results to the Summary Directory..."	
#	cp -v ${REGIONALDIR}/*_nominal.P0_05.txt ${SUMMARY}/
#	cp -v ${REGIONALDIR}/*_perm.P0_05.txt ${SUMMARY}/
#	cp -v ${REGIONALDIR}/*.pdf ${SUMMARY}/
#	
#	echo ""
#	echo "Summarising all results for the ${VARIANT} locus in one file..."
#	
#	echo ""
#	echo "Nominal results..."
#	cat ${SUMMARY}/${TODAY}_ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}_nominal.P0_05.txt | tail -n +2 | awk -v LOCUS_VARIANT=$VARIANT '{ print LOCUS_VARIANT, $0 }' OFS="\t"  >> ${SUMMARY}/ctmm_QC_eqtlnom_summary.txt
#	gzip -v ${SUMMARY}/${TODAY}_ctmm_QC_eqtlnom_${VARIANT}_excl_${EXCLUSION_TYPE}_nominal.P0_05.txt
#	echo ""
#	echo "Permutation results..."
#	cat ${SUMMARY}/${TODAY}_ctmm_QC_eqtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}_perm.P0_05.txt | tail -n +2 | awk -v LOCUS_VARIANT=$VARIANT '{ print LOCUS_VARIANT, $0 }' OFS="\t"  >> ${SUMMARY}/ctmm_QC_eqtlperm_summary.txt
#	gzip -v ${SUMMARY}/${TODAY}_ctmm_QC_eqtlperm_${VARIANT}_excl_${EXCLUSION_TYPE}_perm.P0_05.txt
#	
#done < ${REGIONS}
#
#### ZIPPING FINAL SUMMARY RESULTS
#gzip -v ${SUMMARY}/ctmm_QC_eqtlnom_summary.txt
#gzip -v ${SUMMARY}/ctmm_QC_eqtlperm_summary.txt

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Parse nominal results to get Hits per Locus and (mapped) Gene, and input files"
echo "for LocusZoom v1.2+."

# First we will collect all the nominal association results.
echo ""
echo "Parsing nominal results..."
cd ${SUMMARY}
pwd
#python ${FASTQTLPARSER} ${SUMMARY}/ctmm_QC_eqtlnom_summary.txt.gz

# Now we will start plotting per locus each gene-probe-pair.
echo ""

while IFS='' read -r REGIONOFINTEREST || [[ -n "$REGIONOFINTEREST" ]]; do
	#	1		2		3	4	5		6		7			8		9
	#Variant	Locus	Chr	BP	BP-1Mb	BP+1Mb	WindowSize	Type	Phenotype
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
	echo "        INITIALISING LOCUSZOOM PLOTTING FOR ${VARIANT} LOCUS"
	echo "===================================================================="
	echo ""
	
	# Getting only the top part of the variant-gene-probeid list
	cat ${SUMMARY}/_loci/${VARIANT}.txt | tail -n +2 > ${SUMMARY}/_loci/${VARIANT}.LZ.txt
	LOCUSHITS=${SUMMARY}/_loci/${VARIANT}.LZ.txt
	echo "These are the hits we're interested in for ${VARIANT}..."
	echo ""
	cat ${LOCUSHITS}
	
	echo ""
	while IFS='' read -r VARIANTGENEPROBE || [[ -n "$VARIANTGENEPROBE" ]]; do
		#	1		2			3		4			5
		#	Locus	GeneName	ProbeID	N_Variants	N_Significant
		LINE=${VARIANTGENEPROBE}
		LOCUSVARIANT=$(echo "${LINE}" | awk '{print $1}')
		GENENAME=$(echo "${LINE}" | awk '{print $2}')
		PROBEID=$(echo "${LINE}" | awk '{print $3}')
		N_VARIANTS=$(echo "${LINE}" | awk '{print $4}')
		N_SIGNIFICANT=$(echo "${LINE}" | awk '{print $5}')
	
		echo "===================================================================="
		echo "Plotting results for the ${LOCUSVARIANT} locus on ${CHR}:${START}-${END}."
		echo "	* Plotting association results for ${GENENAME} and ${PROBEID}."
		echo "	* Total number of variants analysed: 	${N_VARIANTS}."
		echo "	* Total number of significant variants:	${N_SIGNIFICANT}."
		
		echo ""
		### FOR DEBUGGING
		#echo "Checking existence of proper input file..."
		#ls -lh ${SUMMARY}/_probes/${LOCUSVARIANT}_${GENENAME}_${PROBEID}.lz
		#echo ""
		#echo "Head:"
		#head ${SUMMARY}/_probes/${LOCUSVARIANT}_${GENENAME}_${PROBEID}.lz
		#echo ""
		#echo "Tail:"
		#tail ${SUMMARY}/_probes/${LOCUSVARIANT}_${GENENAME}_${PROBEID}.lz
		#echo ""
		#echo "Row count:"
		#cat ${SUMMARY}/_probes/${LOCUSVARIANT}_${GENENAME}_${PROBEID}.lz | wc -l
		
		### Setting up LocusZoom v1.2+ plotting
		# Some general settings
		LOCUSZOOM_SETTINGS="ldColors=\"#595A5C,#4C81BF,#1396D8,#C5D220,#F59D10,red,#9A3480\" showRecomb=TRUE dCol='r^2' drawMarkerNames=FALSE refsnpTextSize=0.8 geneFontSize=0.6 showRug=FALSE showAnnot=FALSE showRefsnpAnnot=TRUE showGenes=TRUE clean=TRUE bigDiamond=TRUE ymax=12 rfrows=10 refsnpLineWidth=2"
		# The proper genome-build
		LDMAP="--pop EUR --build hg19 --source 1000G_March2012"
		# Directory prefix
		PREFIX="${LOCUSVARIANT}_${GENE}_${PROBEID}_excl_${EXCLUSION_TYPE}_"
		
		# Actual plotting
		${LZ13} --metal ${SUMMARY}/_probes/${LOCUSVARIANT}_${GENE}_${PROBEID}.lz --markercol MarkerName --pvalcol P-value --delim tab --chr ${CHR} --start ${START} --end ${END} ${LDMAP} ${LOCUSZOOM_SETTINGS} --prefix=${PREFIX} theme=publication title="${LOCUSVARIANT} - ${GENENAME} (${PROBEID})" 
	
	done < ${LOCUSHITS}
	
	rm -v ${SUMMARY}/_loci/${VARIANT}.LZ.txt
	
done < ${REGIONS}

#done

echo ""
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "Wow. I'm all done buddy. What a job! let's have a beer!"
date



