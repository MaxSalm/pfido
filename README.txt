PFIDO  (Phase-Free Inversion Detection Operator)

Version: 1
Freeware 

Author:   Max Salm
Email:    max.salm05@imperial.ac.uk, maxsalm3@gmail.com  


Reference
=========
Genome Res. 2012 Jun;22(6):1144-53. doi: 10.1101/gr.126037.111. Epub 2012 Mar 7.
The origin, global distribution, and functional impact of the human 8p23 inversion polymorphism.
Salm MP, Horswell SD, Hutchison CE, Speedy HE, Yang X, Liang L, Schadt EE, Cookson WO, Wierzbicki AS, Naoumova RP, Shoulders CC.


Contents
========
1. System requirements
2. Tutorial
3. Startup scripts


1. System requirements
======================
PFIDO requires R Version 2.9.2 or higher which is available at http://cran.r-project.org/
The following R packages are also required:

a) snpMatrix	
Clayton, D.G. and Leung, Hin-Tak (2007) An R package for analysis of whole-genome association studies. Human Heredity 64:45-51.

b) mclust
Chris Fraley and Adrian E. Raftery (2006) MCLUST Version 3 for R: Normal Mixture Modeling and Model-based Clustering.Technical Report No. 504, Department of Statistics, University of Washington(revised 2009)
Chris Fraley and Adrian E. Raftery (2002) Model-based Clustering, Discriminant Analysis and Density Estimation. Journal of the  American Statistical Association 97:611-631

c) clValid
Guy Brock, Vasyl Pihur, Susmita Datta and and Somnath Datta (2008). clValid: Validation of Clustering Results. R package version 0.5-7. 

d) extremevalues
Mark van der Loo (2009). extremevalues: Outlier detection in onedimensional data. R package version 1.0.

e) moments
Lukasz Komsta and Frederick Novomestky (2007). moments: Moments, cumulants, skewness, kurtosis and related tests. R package version 0.11. 

f) RColorBrewer
Erich Neuwirth (2007). RColorBrewer: ColorBrewer palettes. R package version 1.0-2.


These will be downloaded and installed upon first use of PFIDO, apart from mclust which has to be downladed manually from:

http://cran.r-project.org/web/packages/mclust/index.html 



2. Instructions
===========
	1) Open and R session
	2) Change Directory to the PFIDO_gui directory (either by selecting File>Change dir.., or typing setwd("<full path>/PFIDO_gui"))
	3) Load PFIDO, by typing source("scripts/PFIDO_gui.R")
	4) Follow GUI instructions
	5) Results are found in the folder "/PFIDO_gui/output"

Either download chr8 HapMap genotypes from ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/latest/forward/non-redundant/and put in the "tmp" folder and unzip, or the routine will connect and download it each time its run.

To check the package runs on your computer, please try to run the example HapMap-derived files (TSI.ped/TSI.map) in "/PFIDO_gui/input"


PFIDO also exists as an R function, and can be loaded via the command: source("PFIDO_module.R")
This allows finer control of the algorithm's parameters than offered by the batch file above.

A. Input files:
	PED/MAP file format (tab-delimmited, see http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped for an example). 
		I recommend using PLINK to correctly prepare the input files. 	
	Family and individual identifiers must be numerically encoded - Individuals should be unrelated and of comparable ancestry.
	Only simple phenotype codes (e.g. 0/1/2/-9) are supported at this time. 
	Alleles are expected to be encoded by A,C,G or T and delimmited by a white-space. 
	Missing genotypes must be encoded by 0. 


B. Options:
	Eighteen options are given.

PED = 8p23_Eur1kgp.ped		<--- Input ped file (required)
MAP = 8p23_Eur1kgp.map		<--- Input map file (required) 
REFERENCE = CEU			<--- Reference HapMap population (CEU/YRI/JPT - required) 
HWE = 0.001			<--- Exclude markers that fail the Hardy-Weinberg test at a specified threshold (p) 
MAF = 0.01			<--- Exclude markers with MAF less than a threshold (0-1: e.g. 0.05) 
call.rate = 0.9			<--- Exclude markers with MAF less than a threshold (0-1: e.g. 0.05) 
mind      = 0.9			<--- Exclude samples missing too much genotype data (0-1: e.g. 0.9 (i.e.90%))
THRESHOLD = 0.05		<--- PFIDO p-value threshold
out.nom   = 8p23_Eur1kgp_output.txt 		<--- PFIDO output filename
restrict  = Y/N			<--- Restrict analysis to optimised SNP sets (see Supplemental Note in paper)
OUTLIER   = Y/N			<--- Identify and remove outliers 

c. Known problems
	1) Please ensure that an appropriate reference population is selected, otherwise the inversion-type calls will likely be completely inaccurate.

