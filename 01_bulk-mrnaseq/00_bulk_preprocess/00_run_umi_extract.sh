#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe smp 2     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N umitools       # Specify job name

conda activate umitools

file_prefix=$1
in_DIR="/afs/crc/group/TIMELab/SpaceU87_Nov2024/ILMN_2346_HarkerBrent_ND_1Pool_Nov2024/"
out_DIR="/afs/crc/group/TIMELab/SpaceU87_Nov2024/UMITOOLS_out/"

umi_tools extract --stdin="${in_DIR}/${file_prefix}_R1_001.fastq.gz" \
		--bc-pattern2=NNNNNNNN \
		--read2-in="${in_DIR}/${file_prefix}_R2_001.fastq.gz" \
		--log="${i}_extract.log" \
		--stdout="${out_DIR}/${file_prefix}_R1_001.fastq.gz" \
		--read2-out="${out_DIR}/${file_prefix}_R2_001.fastq.gz" 
