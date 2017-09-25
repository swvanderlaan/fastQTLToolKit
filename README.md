fastQTLToolKit
============
[![DOI](https://zenodo.org/badge/55601648.svg)](https://zenodo.org/badge/latestdoi/55601648)

This repository contains various scripts in BASH and Python scripts for genetic analyses of (genome-wide) methylation and expression data ('quantitative trait locus', QTL analyses) of the Athero-Express Genomics Studies 1 and 2 (AEGS) or CTMM Genomics Study (CTMM). 

To perform these QTL analyses we make use of [fastQTL](http://fastqtl.sourceforge.net) developed by the lab of [Prof. Emmanouil (Manolis) Dermitzakis](http://funpopgen.unige.ch) and published in [Bioinformatics](https://www.ncbi.nlm.nih.gov/pubmed/26708335). 

AEGS contains methylation data using the Illumina Human Methylation 450K BeadChip from carotid plaques (n = 664). CTMM contains expression data using the Illumina Human HT12 v4r2 BeadChip from monocytes (n = 308). AEGS1 was genotyped using Affymetrix SNP 5.0, AEGS2 using Affymetrix Axiom CEU, and CTMM using Affyemtrix Axiom TX (a custom version of the 'Biobank'-chip). Both AEGS and CTMM were imputed using 1000G phase 3 version 5 and GoNL5 as a reference.

All scripts are annotated for debugging purposes - and future reference. Scripts will work within the context of a certain Linux environment (in this case a CentOS7 system on a SUN Grid Engine background). 

The installation procedure is quite straightforward, and only entails two steps consisting of command one-liners that are *easy* to read. You can copy/paste each example command, per block of code. For some steps you need administrator privileges. Follow the steps in consecutive order.

```
these `mono-type font` illustrate commands illustrate terminal commands. You can copy & paste these.
```

To make it easier to copy and paste, long commands that stretch over multiple lines are structered as follows:

```
Multiline commands end with a dash \
	indent 4 spaces, and continue on the next line. \
	Copy & paste these whole blocks of code.
```

Although we made it easy to just select, copy and paste and run these blocks of code, it is not a good practise to blindly copy and paste commands. Try to be aware about what you are doing. And never, never run `sudo` commands without a good reason to do so. 

We have tested fastQTLToolKit on CentOS7, OS X El Capitan (version 10.11.[x]), and macOS Sierra (version 10.12.[x]). 

--------------

#### Requirements

* ***R Statistics*** -- The scripts are tested in conjunction with 3.3.1+. You can find R here: https://www.r-project.org.
* ***For R*** -- These libraries and their dependencies are needed for the R-script to work. 
	- "optparse"
	- "tools"
	- "qvalue"
* ***QCTOOL*** -- Needed for the extraction and conversion of required genotype data. Required: v1.5. Download: http://www.well.ox.ac.uk/~gav/resources/qctool_v1.5-linux-x86_64.tgz.
* ***SNPTEST*** -- Needed for the calculation of summary statistics of the extracted genotype data. Required: v2.5.2. Website: https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html.
* ***fastQTL*** -- Needed to run fastQTL. Required: v2.184. Website: http://fastqtl.sourceforge.net.
* ***LocusZoom*** -- Needed to generate LocusZoom plots. Required: v1.3. Website: http://genome.sph.umich.edu/wiki/LocusZoom_Standalone.
* ***BGZip*** -- Needed for VCF conversion and gzipping. Required: as part of the htslib-1.3 package. Website: http://www.htslib.org/download/.
* ***Tabix*** -- Needed for VCF indexing. Required: as part of the htslib-1.3 package. Website: see link above.

####_NOTE: the above scripts, tools, and software are all installed on our local system._

--------------

#### Installing the scripts locally

You can use the scripts locally to run analyses on a Unix-based system, like Mac OS X (Mountain Lion+). We need to make an appropriate directory to download 'gits' to, and install this 'git'.

##### Step 1: make a directory, and go there.

```
mkdir -p ~/git/ && cd ~/git
```

##### Step 2: clone this git, unless it already exists.

```
if [ -d ~/git/GWASToolKit/.git ]; then \
		cd ~/git/GWASToolKit && git pull; \
	else \
		cd ~/git/ && git clone https://github.com/swvanderlaan/fastQTLToolKit.git; \
	fi
```


--------------

#### Quantitative Trait Locus analyses 

You can select the type of analysis by providing the following _obligatory_ 11 arguments. Some relevant statistics, such as HWE, minor allele count (MAC), and coded allele frequency (CAF) will be added to the final summarized result. LocusZoom style figures will be made automatically for _eQTL-analyses_ alone. 


* ***Argument #1*** -- indicate which study type you want to analyze, so either [AEMS450K1/AEMS450K2/CTMM]:
	- ***AEMS450K1*** -- methylation quantitative trait locus (mQTL) analysis on plaques (n = 444) or blood (n = 93) in the Athero-Express Methylation Study 450K phase 1.
	- ***AEMS450K2*** -- mQTL analysis on plaques in the Athero-Express Methylation Study 450K phase 2 (n = 158).
	- ***CTMM*** -- expression QTL (eQTL) analysis in monocytes from CTMM.
* ***Argument #2*** -- the sample type must be [AEMS450K1: PLAQUES/BLOOD], [AEMS450K2: PLAQUES], or [CTMM: MONOCYTES].
* ***Argument #3*** -- the root directory, e.g. /server/folderX/folderY/folderZ/someqtlanalysis.
* ***Argument #4*** -- where you want stuff to be save inside the rootdir,  e.g. mqtl_aems450k1
* ***Argument #5*** -- project name, e.g. 'CAD'.
* ***Argument #6*** -- text file with on each line the regions of interest, refer to example files.
* ***Argument #7*** -- the type of exclusion to apply:
	- ***AEMS/CTMM*** -- DEFAULT/SMOKER/NONSMOKER/MALES/FEMALES/T2D/NONT2D
	- ***AEMS-specific*** -- CKD/NONCKD/PRE2007/POST2007/NONAEGS/NONAEGSFEMALES/NONAEGSMALES.
* ***Argument #8*** -- text file with excluded covariates, refer to example file.
* ***Argument #9*** -- qsub e-mail address, e.g. yourname@mailadres.com.
* ***Argument #10*** -- qsub mail settings, e.g. 'beas' - refer to qsub manual.
* ***Argument #11*** -- configuration file.

--------------

#### Exclusion options

####_CTMM/AEGS_

- ***DEFAULT*** -- include all CTMM or AEGS samples (note: _only_ carotid endartectomy [CEA] patients are considered for AEGS).
- ***SMOKER*** -- include only individuals with 'current smoking' status = 'yes'.
- ***NONSMOKER*** -- include only individuals with 'current smoking' status = 'no'.
- ***MALES*** -- include only 'female' individuals (exclude 'males').
- ***FEMALES*** -- include only 'male' individuals (exclude 'females').
- ***T2D*** -- include only individuals _with_ 'type 2 diabetes' based on diagnosis and/or medication use.
- ***NONT2D*** -- include only individuals _without_ 'type 2 diabetes' based on diagnosis and/or medication use.

####_AEGS specific_
- ***CKD*** -- include only individuals _with_ 'chronic kidney disease' (CKD>2) based on eGFR (using creatinine) using KDOQI criteria.
- ***NONCKD*** -- include only individuals _without_ 'chronic kidney disease' (CKD>2) based on eGFR (using creatinine) using KDOQI criteria.
- ***PRE2007*** -- include only individuals operated _before_ 2007.
- ***POST2007*** -- include only individuals operated _after_ 2007.
- ***NONAEGS*** -- include all also non-CEA individuals.
- ***NONAEGSMALES*** -- include all also non-CEA individuals, but only 'female' individuals (exclude 'males').
- ***NONAEGSFEMALES*** -- include all also non-CEA individuals, but only 'male' individuals (exclude 'females').

--------------

#### File descriptions

- ***excl_cov_eqtl.txt*** -- list of covariates that should be excluded from analysis (for CTMM analysis).
- ***excl_cov_sex_eqtl.txt*** -- list of covariates that should be excluded from analysis, when running sex-stratified analyses (for CTMM analysis).
- ***excl_cov.txt*** -- list of covariates that should be excluded from analysis (for AEGS analysis).
- ***excl_cov_sex.txt*** -- list of covariates that should be excluded from analysis, when running sex-stratified analyses (for AEGS analysis).
- ***regions_for_qtl.small.txt*** -- example list of regions (short list).
- ***regions_for_qtl.txt*** -- example list of regions.
- ***fastQTLAnalyzer.sh*** -- Main analysis-script.
- ***fastQTLChecker.sh*** -- Script that checks the analysis output.
- ***fastQTLSummarizer.sh*** -- Script that summarizes the fastQTL analysis results.
- ***fastQTL_QC.R*** -- R script that does some parsing of results and QC, as well as generating some QC plots.
- ***NominalResultsParser.py*** -- Script need to parse the nominal results data to obtain per-locus lists of number of variants, transcripts, and the LocusZoom input file.
- ***fastQTLPlotter.sh*** -- Script that generates per transcript (!) for each gene and locus a LocusZoom style plot (CTMM analysis only!).
- ***parse_clumps_eqtl.pl*** -- Clumps results per locus -- NOTE: BETA-version; not used yet.

--------------

#### Roadmap
fastQTLToolKit will no longer be updated, rather a new QTLToolKit is underway which is based on [QTLTool](https://www.nature.com/articles/ncomms15452) available [here](https://qtltools.github.io/qtltools/). Below a list of features/issues which will handled in QTLToolKit.
- update NominalResultsParser for mQTL analyses.
- generate regional plots for mQTL analyses.
- generalisation: modify 'exclusion type'-part to work with any exclusion list for a particular dataset.
- generalisation: modify 'input'-data to work with any particular dataset.
- add in the existing 'make-bed-files' scripts.
- add in trans-QTL analyses pipeline.

--------------

#### The MIT License (MIT)
####Copyright (c) 2015-2017 Sander W. van der Laan | s.w.vanderlaan-2 [at] umcutrecht.nl

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:   

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.
