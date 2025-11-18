#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe smp 8     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N ftrCounts       # Specify job name

file_prefix=$1
conda activate featurecounts
in_DIR="/afs/crc/group/TIMELab/SpaceU87_Nov2024/STAR_out"

featureCounts -s 2 \
      -p \
      -a /afs/crc/group/TIMELab/genomes/hsapiens/refseq/ncbi_dataset/data/GCF_000001405.40/genomic.gtf \
      -o /afs/crc/group/TIMELab/SpaceU87_Nov2024/featureCounts/${file_prefix}_counts.txt \
      "${in_DIR}/${file_prefix}/${file_prefix}_deduplicated.bam"
