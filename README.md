# IGOFingerPrinting: Tumor/Normal pair analysis using IGOFingerPrinting commandline tool

This software uses BAM files (without duplicates marked) as input. The libraries for samples are prepared using a SNP panel to sequence the target regions around 39 SNPs used for the analysis.  The software scripts and resources are located under /ifs/res/GCL/hiseq/Stats/IGOFingerPrinting/ . You can copy this directory to any directory of your choice and replace /ifs/res/GCL/hiseq/Stats/ with path to new directory in the examples below. To run the analysis you will need the following

1.	Bed File containing coordinates for target SNPS. A bed file is already created under  target_snps directory inside the software directory. This file needs to be created only once every time list of target SNPs change.
2.	VCF File containing detailed information about the SNPS. This file will be used by GATK software to insert SNP ID's to variant call files generated after the analysis. This file needs to be updated only once every time target SNPs change.
3.	Patient IDs File Under vialelab/LIMS/FingerPrinting/  create a directory named as PROJECT_ID that need analysis. If it is single project then it should be named like 08494_B. If you need multiple projects to be analyzed together, then name the directory with PROJECT_ID's separated by "_". 
example: 06000_EA_08822_E_08822_F. 
Inside the directory, create a tab delimited text file with the same name as parent directory and add IGO ID, CMO ID, PATIENT ID separated by "tab" for each sample. If needed, you can refer to some of the files I have already created. Finally, save this text file as "Unix LF" format. 
Note: It is important that the directory name and the file name are same except the file extension.

# Analysis

To run the analysis run the following command from /ifs/res/GCL/hiseq/Stats/IGOFingerPrinting/

# Example:
python RunFingerPrintingAnalysis.py -run_project PITT_0251_AHW5KYBBXX--06000_EA,PITT_0251_AHW5KYBBXX--08822_E,PITT_0251_AHW5KYBBXX--08822_F -replace_ids 06000_EA_08822_E_08822_F
You need to provide two values to the script

-run_project 	For each project provide the RUNID and ProjectID separated by "--". If there are multiple projects for analysis, or a project was sequenced over multiple runs, separate the RUNID-PROJECTID by ",".
	Example:  PITT_0251_AHW5KYBBXX--06000_EA
or
	PITT_0251_AHW5KYBBXX--06000_EA,PITT_0251_AHW5KYBBXX--08822_E
	
Note: It is important that the ProjectIds are always prefixed with Zero.

-replace_ids 	A Patient ID file containing Patient Ids for all the samples being analyzed. This file should be stored as described above under step 3.


# Results:

Under the Software directory, there is a directory called "reports". The results are stored under this directory. The software automatically creates the result directory. The software generates following outputs that should be shared as results

"	Heatmap pdf: file with sample pairs if found. 
Eg: 05732_AA_05732_AB_05732_X_heatmap.pdf
"	Phylogenetic tree pdf : file with sample pairs if found.
Eg: 05732_AA_05732_AB_05732_X_tree.pdf
"	P-values: file containing match and pvalues annotation for samples.
Eg: 05732_AA_05732_AB_05732_X_p-values.csv
"	VCF: Variant Call File containing the Genotype info for the target SNP's.
Eg: 06000_EA_08822_E_08822_F_with_patientid.vcf

# Intermediate Results files:

Some intermediate files are also created during analysis that could be useful for troubleshooting. 
"	VCF: file originally created during the analysis containing Sample Ids extracted from bam files. These samples are then changed to corresponding patient ID's using the Patient ID file provided in the command.
"	Genotype Numeric Matrix: For correlation analysis the Genotypes for each SNP for each sample are converted to number and output is stored in a file. Eg: genotype_num_matrix.txt 
"	Bam input list: file containing path to the bam files for all the sampes being analyzed. Eg: bam_input.list

# Additional Details:

The software utilizes HaplotyperCaller from GATK (version 3.7) software for making the variant calls. The variant calls are then converted into a numeric metrix using following dictionary:

NT_TO_NUMBER = {'A': '1', 'T': '2', 'G': '3', 'C': '4', '.': '0'}

For example: a genotype 'A/G' in the Variant file is converted to 14, 'C/G' is converted to 43, ./. is converted to 00 and G/. is converted to 30. 

The numeric matrix is then analyzed in R. Correlation matrix is created with "pearson" method and is used to generate Heatmap, tree, and p-values files.








