#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe smp 4     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N trimmomatic       # Specify job name

conda activate umitools

file_name=$1
in_DIR="/afs/crc/group/TIMELab/SpaceU87_Nov2024/UMITOOLS_out/"
out_DIR="/afs/crc/group/TIMELab/SpaceU87_Nov2024/UMITOOLS_out/"

trimmomatic SE -threads 4 "${in_DIR}/${file_name}_R2_001.fastq.gz" "${out_DIR}/trimmed/${file_name}_R2_001.fastq.gz" HEADCROP:6
