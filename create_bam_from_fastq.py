import glob
import os
import re

def bam_from_fastq(fastq_dir, BWA, SAMTOOLS, REF_HG19):
  fastq_files_r1 = glob.glob(fastq_dir + '/*R1*.fastq.gz')
  for r1 in fastq_files_r1:
    r2 = r1.replace('R1', 'R2')
    out = r1.replace(r1[r1.find("R1"):], "PE.sam")
    sample_name = r1.split("/")[-1].split("_R1_")[0]
    read_group_sample_name = sample_name[::-1][sample_name[::-1].find("S_")+2:][::-1]
    read_info = os.popen('zcat < ' + r1 + ' | head -1').read()
    read_group= re.split(":",read_info)
    machine_run_id = read_group[0][1:] + "_" + read_group[1]
    platform = 'illumina'
    company = "GCL@MSKCC"
    read_group_string = "@RG\\tID:" + machine_run_id + "\\tLB:" + read_group_sample_name + "\\tSM:" + read_group_sample_name + "\\tPL:illumina\\tCN:" + company
    os.system(BWA + ' mem -P -R "%s" %s %s %s > %s' %(read_group_string,REF_HG19,r1,r2,out))

  sam_files = glob.glob(fastq_dir + '/*.sam')
  for sam in sam_files:
    bam_out = sam.replace('.sam','.bam')
    os.system(SAMTOOLS + ' view -bT %s %s > %s' %(REF_HG19,sam,bam_out))
    os.system('rm '+ sam)
    os.system(SAMTOOLS + ' index -b %s' %(bam_out))
  
  bam_files = glob.glob(fastq_dir + '/*.bam')
  for bam in bam_files:
    sorted_bam = bam.replace("PE.bam","PE_sorted.bam")
    temp_sorted_bam = bam.replace(".bam","_aln_sorted.bam")
    os.system(SAMTOOLS + ' sort -T %s -o %s %s' %(temp_sorted_bam,sorted_bam,bam))
    os.system(SAMTOOLS + ' index -b %s' %(sorted_bam))
    os.system('rm ' + bam)
