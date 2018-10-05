#!/opt/common/CentOS_6/python/python-2.7.8/bin/
# coding: utf-8

import argparse
import glob
import logging
import os
import re
import sys
import csv
import vcf

logging.basicConfig(level=logging.INFO,
                    format='%(levelname)-5s @ %(asctime)s:\n\t %(message)s \n',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )
error = logging.critical
warn = logging.warning
debug = logging.debug
info = logging.info

__version__ = "0.0.1"

parser = argparse.ArgumentParser(prog="IGOSampleMatchMaker",
                                 usage="A tool to match tumor/normal pairs using list of SNPs")
parser.add_argument('-run_project', type=str,
                    help='Provide a comma separated list of run-project to analyze. Eg: '
                         'A00227_0009_AH2L7YDMXX-06000, DIANA_0101_A12KJDMXX-06000_BS.Run ID and Project number '
                         'must be separated by "-".',
                    required=True)
parser.add_argument('-replace_ids', type=str,
                    help='If need to replace sample Ids in vcf file with new Ids containing patient Id, use this '
                         'option. Patient Id file should be a tab delimited without header with data in order '
                         'sampleID, cmoID and patientID.',
                    required=False)

args = parser.parse_args()
projects = str(args.run_project)
patient_id_file_name = str(args.replace_ids)
matchmaker_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
JAVA = "/opt/common/CentOS_6/java/jdk1.8.0_31/bin/java"
PICARDS = JAVA + " -jar /opt/common/CentOS_6/picard/picard-2.9.0-jar/picard.jar"
GATK = JAVA + " -jar /opt/common/CentOS_6/gatk/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar"
BWA = '/opt/common/CentOS_6/bwa/bwa-0.7.12/bwa'
SAMTOOLS = '/opt/common/CentOS_6/samtools/samtools-1.2/samtools'
# GATK = JAVA + " -jar " + matchmaker_dir + "/gatk-4.0.2.1/gatk-package-4.0.2.1-local.jar"
R = "/opt/common/CentOS_6/R/R-3.2.0/bin/Rscript"
REF = "/ifs/data/bio/Genomes/H.sapiens/hg19/human_hg19.fa"
REF_BWA = "/ifs/data/bio/Genomes/H.sapiens/hg19/BWA_0.7.5a/human_hg19.fa"
TARGET_SNP_FILE = matchmaker_dir + "/target_snps/Fingerprinting_39_snps.bed"
file_with_patient_ids = matchmaker_dir + "/Sample_ID_Files/" + patient_id_file_name
SNP_VCF = matchmaker_dir + "/target_snps/dbsnplist_snp.vcf"
#FASTQ_INPUT_DIR_SERVER = "/ifs/input/GCL/hiseq/FASTQ/"
INPUT_DIR_SERVER = "/ifs/res/GCL/hiseq/Stats/"
NT_TO_NUMBER = {'A': '1', 'T': '2', 'G': '3', 'C': '4', '.': '0'}

if not projects:
    raise Exception("missing arguement -run_project")


#def copy_fastq_files(project_list, server_fastq_dir, destination_directory):
#    run_project = project_list.split(',')
#    for item in run_project:
#        run_project_info = item.split('-')
#        run_id = run_project_info[0].strip()
#        project = 'Project_' + run_project_info[1].strip()
#        project_fastq_source_dir = server_fastq_dir + run_id + '/' + project
#        os.system('cp ' + project_fastq_source_dir + '/Sample_*/*.fastq.gz ' + destination_directory)
#        os.system("echo Done copying FASTQ files for " + project)


#def create_bam_input_file(bam_file_dir, reports_dir):
#    input_bam_list_file = reports_dir + "/bam_input.list"
#    if os.path.exists(input_bam_list_file):
#        os.remove(input_bam_list_file)

#    project_bam_file_names = glob.glob(bam_file_dir + "/*sorted.bam")
#    with open(input_bam_list_file, "a") as bam_file:
#        for filename in project_bam_file_names:
#            bam_file.write(os.path.abspath(filename + "\n"))
#    bam_file.close()
#    os.chmod(input_bam_list_file, 0o777)
#    os.system("echo Done Creating file with path to bam files.")
#    return input_bam_list_file
#


def create_bam_input_file(projects_list, reports_dir):
    input_bam_list_file = reports_dir + "/bam_input.list"
    if os.path.exists(input_bam_list_file):
        os.remove(input_bam_list_file)
    projects_list = projects_list.split(',')
    file_bam_list = open(input_bam_list_file, "a")
    for item in projects_list:
        os.system("echo " + item)
        run_id = get_run_id(item)
        project_id = get_project_id(item)
        project_bam_file_names = glob.glob(INPUT_DIR_SERVER + run_id + "/*___P" + project_id + "___*_FP_*.bam")
        for filename in project_bam_file_names:
            print  filename
            file_bam_list.write(os.path.abspath(filename + "\n"))
    file_bam_list.close()
    os.system("chmod 0775 " + input_bam_list_file)
    os.system("echo Done Creating file with path to bam files.")
    return input_bam_list_file


def get_reports_dir_name(project_list):
    run_project = project_list.split(',')
    dir_name = []
    if len(run_project) == 1:
        return run_project[0].split('-')[1].strip()
    if len(run_project) < 1:
        raise Exception(
            "Cannot create Project directory in fastq_files directory. Reason: Cannot find project ID's in provided "
            "run-project inputs.")
    for item in run_project:
        project_id = item.split('-')[1].strip()
        if project_id not in dir_name:
            dir_name.append(project_id)
    return '_'.join(dir_name)


def get_run_id(project):
    return project.split('-')[0]


def get_project_id(project):
    return project.split('-')[1]


def create_project_directory(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
        os.system("chmod 0775 " + directory_name)


def get_lines_with_genotypes_from_vcf(vcf_file):
    genotype_info = []
    if os.stat(vcf_file).st_size == 0:
        raise Exception("VCF file " + vcf_file + " is empty.")
    with open(vcf_file, "r") as infile:
        for line in infile:
            if "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT" in line or line.startswith(
                    "chr"):
                genotype_info.append(line)
    if len(genotype_info) == 0:
        raise Exception("Provided VCF file " + vcf_file + "  does not contain any genotype data.")
    return genotype_info


def allele_to_number(allele, ref, alt):
    allele1 = None
    if allele == '0':
        return str(NT_TO_NUMBER[ref])
    if allele == '1':
        return str(NT_TO_NUMBER[alt])
    if allele == '.':
        return str(NT_TO_NUMBER['.'])
    if allele != '.\n' and allele != '.' and int(allele) > 1:
        return str(NT_TO_NUMBER[alt.split("/")[allele - 1]])
    else:
        return str(NT_TO_NUMBER['.'])


def change_sample_ids(vcf_header, patientid_file):
    with open(patientid_file, 'r') as patient_data:
        for line in patient_data:
            sample_data = line.split("\t")
            sample_id = sample_data[0]
            cmo_id = sample_data[1]
            patient_id = sample_data[2]
            new_id = sample_id.strip('\n') + "-" + patient_id.strip('\n')
            for header_value in vcf_header:
                if sample_id in header_value and cmo_id in header_value:
                    vcf_header[vcf_header.index(header_value)] = new_id
    patient_data.close()
    return "\t".join(vcf_header) + "\n"


def replace_sampleid_with_patientid(vcf_file, patientid_file):
    vcf_with_patientid_in_header = str(os.path.dirname(os.path.abspath(vcf_file))) + "/" + str(
        os.path.splitext(os.path.basename(vcf_file))[0]) + "_with_patientid.vcf"
    if os.stat(vcf_file).st_size == 0:
        raise Exception("VCF file " + vcf_file + " is empty.")
    if os.stat(patientid_file).st_size == 0:
        raise Exception("Patient ID file " + patientid_file + " is empty.")
    if os.path.exists(vcf_with_patientid_in_header):
        os.remove(vcf_with_patientid_in_header)
    with open(vcf_file, 'r') as vcf, open(vcf_with_patientid_in_header, "a") as new_vcf:
        for line in vcf:
            if line.startswith("#CHR"):
                vcf_header_values = line.split("\t")
                vcf_header_with_patient_ids = change_sample_ids(vcf_header_values, patientid_file)
                new_vcf.write(vcf_header_with_patient_ids)
            else:
                new_vcf.write(line)
    vcf.close()
    new_vcf.close()
    os.system("chmod 0775 " + vcf_with_patientid_in_header)
    return vcf_with_patientid_in_header


def display_fewer_reads_warnings(variant_data, snp_id, header, position):
    snp_data = re.split(':', variant_data)
    if len(snp_data) == 1:
        #print "Fingerprinting WARNING: snp calls not available for SNP '%s' for sample '%s'" % (
        #    snp_id, header[position])
        os.system("echo Fingerprinting WARNING: snp calls not available for SNP " + snp_id + " for sample " + header[position])
    if len(snp_data) > 1:
        num_reads_ref_allele, num_reads_alt_allele = snp_data[1].split(',')
        if int(num_reads_ref_allele) < 10 and int(num_reads_alt_allele) < 10:
            #print  "Fingerprinting WARNING: Less than 10 reads per allele found for SNP '%s' for sample '%s'" % (
               # snp_id, header[position])
            os.system("echo Fingerprinting WARNING: Less than 10 reads per allele found for SNP " + snp_id + " for sample " + header[position])  

def vcf_to_numeric_matrix(vcf_file, output_directory):
    matrix_file = output_directory + "/genotype_num_matrix.txt"
    if os.path.exists(matrix_file):
        os.remove(matrix_file)

    data_lines = get_lines_with_genotypes_from_vcf(vcf_file)
    header = data_lines[0].split("\t")[9:]
    with open(matrix_file, "a") as matrix:
        matrix.write("\t".join(header) + "\n")
        for line in data_lines:
            if line.startswith("chr"):
                data = line.split("\t")
                snp_data = data[9:]
                snp_id = data[2]
                ref = data[3]
                alt = data[4]
                snp_data_for_all_samples = snp_id
                for sample_data in snp_data:
                    position = snp_data.index(sample_data)
                    display_fewer_reads_warnings(sample_data, snp_id, header, position)
                    genotype = re.split(':', sample_data)[0]
                    allele1, allele2 = genotype.split('/')
                    allele1_to_number = allele_to_number(allele1, ref, alt)
                    allele2_to_number = allele_to_number(allele2, ref, alt)
                    snp_data_for_all_samples = snp_data_for_all_samples + "\t" + allele1_to_number + allele2_to_number
                matrix.write(snp_data_for_all_samples + "\n")
    matrix.close()
    os.system("chmod 0775 " + matrix_file)
    return matrix_file


def annotate_results_for_match(csvfile):
    file1 = open(csvfile, 'rb')
    reader = csv.reader(file1)
    new_rows_list = []
    new_header = ['', 'row', 'column', 'cor', 'p', 'sample1','sample2','match']
    new_rows_list.append(new_header)
    for row in reader:
        if row[1]=="row":
            continue
        new_row = [row[0], row[1], row[2], row[3], row[4], row[1].split('_')[-1], row[2].split('_')[-1]]
        if new_row[5] == new_row[6] and float(new_row[3]) > 0.8 and float(new_row[4]) < 0.05:
           new_row.append('OK')
        else:
           new_row.append('FALSE')
        new_rows_list.append(new_row)
    file1.close()

    with open(csvfile, "wb") as file:
        writer = csv.writer(file)
        writer.writerows(new_rows_list)
    file.close()

def get_genotypes_from_vcf(vcf_file):
    VCF_FILE = vcf_file
    SNP_FP = list() # initialize list
    VCF = vcf.Reader(open(VCF_FILE, 'r'))  # open VCF file
    for RECORD in VCF:
        RECORD_DATA = list()
        RECORD_DATA = [RECORD.CHROM,RECORD.POS,RECORD.REF,RECORD.ALT]
        for SAMPLE in RECORD.samples:
            RECORD_DATA.append(SAMPLE.gt_bases)
        SNP_FP.append(RECORD_DATA)

    SAMPLE_HEADS = list() # create sample headings
    for SAMPLE in RECORD.samples:
        SAMPLE_HEADS.append(SAMPLE.sample)
        
    RECORD_HEAD = ['CHROM','POS','REF','ALT'] # add additional headings to the sample headings
    RECORD_HEAD.extend(SAMPLE_HEADS) # combine record and sample headings
    SNP_FP.insert(0, RECORD_HEAD) # insert headings into the SNP_FP list in the 0 position
    #print(SNP_FP) # sanity check
    CSV_STR = VCF_FILE[:-4] + '_genotypes.csv' # create ouput file name          
    CSV_FILE = open(CSV_STR,'w') # open output file    

    with CSV_FILE: # write SNP_FP list to the output file
        writer = csv.writer(CSV_FILE)
        writer.writerows(SNP_FP)
    os.system("echo Done creating Genotypes file.") 

def main():
    FINGERPRINTING_REPORT_DIR = matchmaker_dir + "/reports/"
    REPORTS_DIR = get_reports_dir_name(projects)
    #INPUT_FASTQ_BAM_DIR = matchmaker_dir + "/input_files/" + REPORTS_DIR
    PROJECT_REPORT_DIR = FINGERPRINTING_REPORT_DIR + "Project_" + REPORTS_DIR

    # Create project directory for reports
    create_project_directory(PROJECT_REPORT_DIR)
    os.system("echo " +  PROJECT_REPORT_DIR)
    # Create a text file with list of all the bam files to be analyzed
    BAM_INPUT_LIST_FILE = create_bam_input_file(projects, PROJECT_REPORT_DIR)
    # BAM_INPUT_LIST_FILE = create_bam_input_file(INPUT_FASTQ_BAM_DIR, PROJECT_REPORT_DIR)
    # BAM_INPUT_LIST_FILE = "/home/sharma/IGOSampleMatchMaker/reports/Project_05732_AA_05732_AB_05732_X/bam_input.list"

    # File to hold vcf output results
    VCF_RESULTS_FILE_NAME = REPORTS_DIR + ".vcf"
    VCF_RESULT_OUTPUT_FILE = PROJECT_REPORT_DIR + "/" + VCF_RESULTS_FILE_NAME

    # Run GATK HaplotypeCaller on the bam files to generate SNP genotypes
    os.system(
        GATK + " -T HaplotypeCaller -stand_call_conf 30 -R " + REF + " -L " + TARGET_SNP_FILE + " -I " + BAM_INPUT_LIST_FILE + " -o " + VCF_RESULT_OUTPUT_FILE + " -mbq 30 -D " + SNP_VCF + ' >> FingerPrinting.log')

    # Add patient ID's to Sample ID's in the vcf file
    PATH_TO_PATIENT_IDS_FILE = matchmaker_dir + "/Sample_ID_Files/" + patient_id_file_name + "/" + patient_id_file_name + ".txt"
    VCF_WITH_PATIENT_ID = replace_sampleid_with_patientid(VCF_RESULT_OUTPUT_FILE, PATH_TO_PATIENT_IDS_FILE)

    # Convert SNP genotypes from vcf file to numeric matrix, which will be used by R script for correlation analysis.
    matrix = None
    matrix = vcf_to_numeric_matrix(VCF_WITH_PATIENT_ID, PROJECT_REPORT_DIR)

    os.system(R + " " + matchmaker_dir + "/phyloTreeGenerator.R --input " + matrix + " --out " + PROJECT_REPORT_DIR +
              " --name " + REPORTS_DIR + ' >> FingerPrinting.log')
    
    # Add match or not match annotation to file with p-values
    P_VALUE_FILE = PROJECT_REPORT_DIR + "/" + REPORTS_DIR + "_p-values.csv"
    annotate_results_for_match(P_VALUE_FILE)
    
    # Extract genotypes from vcf file
    get_genotypes_from_vcf(VCF_RESULT_OUTPUT_FILE)
'''
#Create project directory where project fastq files will be copied and then converted to bam files for input to GATK for fingerprinting analysis. 
    create_project_directory(INPUT_FASTQ_BAM_DIR)
    print INPUT_FASTQ_BAM_DIR

#Copy project bam files to local fastq files directory
    copy_fastq_files(projects, FASTQ_INPUT_DIR_SERVER, INPUT_FASTQ_BAM_DIR)

#Create bam files from copied project fastq files.
    bam_from_fastq(INPUT_FASTQ_BAM_DIR, BWA, SAMTOOLS, REF_BWA)
'''

if __name__ == "__main__":
    main()
