#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe smp 8     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N dedup       # Specify job name

conda activate umitools

file_prefix=$1
in_DIR="/afs/crc/group/TIMELab/SpaceU87_Nov2024/STAR_out/"

cd "${in_DIR}/${file_prefix}"
umi_tools dedup -I "${file_prefix}Aligned.sortedByCoord.out.bam" --output-stats=deduplicated -S "${file_prefix}_deduplicated.bam"
