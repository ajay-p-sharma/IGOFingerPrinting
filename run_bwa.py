import glob
import os
import re

fastq_files_r1 = glob.glob('/home/sharma/NGSCheckMate/fastq_files/Project_05732_AA_05732_AB_05732_X/*R1*.fastq.gz')

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
  os.system('/opt/common/CentOS_6/bwa/bwa-0.7.12/bwa mem -P -R "%s" /ifs/data/bio/Genomes/H.sapiens/hg19/BWA_0.7.5a/human_hg19.fa %s %s > %s' %(read_group_string,r1,r2,out))

sam_files = glob.glob('/home/sharma/NGSCheckMate/fastq_files/Project_05732_AA_05732_AB_05732_X/*.sam')
for sam in sam_files:
  bam_out = sam.replace('.sam','.bam')
  os.system('/opt/common/CentOS_6/samtools/samtools-1.2/samtools view -bT /ifs/data/bio/Genomes/H.sapiens/hg19/BWA_0.7.5a/human_hg19.fa %s > %s' %(sam,bam_out))
  os.system('rm '+ sam)
  os.system('/opt/common/CentOS_6/samtools/samtools-1.2/samtools index -b %s' %(bam_out))
bam_files = glob.glob('/home/sharma/NGSCheckMate/fastq_files/Project_05732_AA_05732_AB_05732_X/*.bam')

for bam in bam_files:
  sorted_bam = bam.replace("PE.bam","PE_sorted.bam")
  temp_sorted_bam = bam.replace(".bam","_aln.sorted")
  os.system('/opt/common/CentOS_6/samtools/samtools-1.2/samtools sort -T /temp/%s -o %s %s' %(temp_sorted_bam,sorted_bam,bam))
  os.system('/opt/common/CentOS_6/samtools/samtools-1.2/samtools index -b %s' %(sorted_bam))
  os.system('rm ' + bam)
