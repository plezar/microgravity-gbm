#!/bin/bash
#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m ae            # Send mail when job begins, ends and aborts
#$ -pe smp 1     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N index       # Specify job name

sample_prefix=$1
DIR="/afs/crc/group/TIMELab/SpaceU87_Nov2024/STAR_out/"
cd "${DIR}/${sample_prefix}"
~/Private/soft/samtools-1.20/samtools index "${sample_prefix}Aligned.sortedByCoord.out.bam"
