#!/opt/ciommon/CentOS_6/python/python-2.7.8/bin/
# coding: utf-8

import argparse
import logging
import os
import sys

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
                                 usage="A tool to match tumor/normal pairs using list of known target SNPs.")
parser.add_argument('-run_project', type=str,
                    help='Provide a comma separated list of run-project to analyze. Eg: '
                         'A00227_0009_AH2L7YDMXX-06000, DIANA_0101_A12KJDMXX-06000_BS.Run ID and Project number '
                         'must be separated by "-".',
                    required=True)
parser.add_argument('-replace_ids', type=str,
                    help='If sample ID in vcf file should contain Patient ID, use this option.'
                         'Provide name without extension of tab delimited without header with data in order '
                         'sampleID, cmoID and patientID.',
                    required=True)

args = parser.parse_args()
run_project_ids = str(args.run_project)
replace_ids_file = str(args.replace_ids)

matchmaker_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
NEW_SAMPLE_ID_FILE_SOURCE = "/pskis34/LIMS/FingerPrinting/" + replace_ids_file + "/" + replace_ids_file + ".txt"
SAMPLE_ID_FILE_DESTINATION = matchmaker_dir + "/Sample_ID_Files/" + replace_ids_file

def main():
    if not os.path.exists(NEW_SAMPLE_ID_FILE_SOURCE):
        os.system("echo Cannot find file with patient Id information. Please make sure that the file is present at " \
                  "/pskis34/LIMS/FingerPrinting/" + replace_ids_file)
        quit()
    if not os.path.exists(SAMPLE_ID_FILE_DESTINATION):
        os.system("mkdir " + SAMPLE_ID_FILE_DESTINATION)
    os.system("echo cp " + NEW_SAMPLE_ID_FILE_SOURCE + " " + SAMPLE_ID_FILE_DESTINATION)
    if not os.path.exists(SAMPLE_ID_FILE_DESTINATION + "/" + replace_ids_file + ".txt"):
        os.system("cp " + NEW_SAMPLE_ID_FILE_SOURCE + " " + SAMPLE_ID_FILE_DESTINATION)
        os.system("chmod 0775 " + SAMPLE_ID_FILE_DESTINATION + "/" + replace_ids_file + ".txt")
    os.system("qsub -N FingerPrinting -cwd -b y -j y /opt/common/CentOS_6/python/python-2.7.8/bin/python " + matchmaker_dir + "/igoMatchMaker.py" + " -run_project " + run_project_ids + " -replace_ids " + replace_ids_file)
    #os.system("qsub -N FingerPrinting -cwd -b y -j y /opt/common/CentOS_6/python/python-2.7.8/bin/python igoMatchMaker.py -run_project " + project_ids + " -replace_ids " + replace_ids)

if __name__ == "__main__":
    main()
