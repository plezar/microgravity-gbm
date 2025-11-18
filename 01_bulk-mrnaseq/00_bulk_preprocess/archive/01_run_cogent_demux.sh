#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe smp 24     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N cogentAP       # Specify job name

in_DIR="/afs/crc/group/TIMELab/SpaceU87_Nov2024/ILMN_2346_HarkerBrent_ND_1Pool_Nov2024/"

$COGENT_AP_HOME/cogent rna demux \
                        -f  "${in_DIR}/combo_R1.fastq.gz" \
                        -p "${in_DIR}/combo_R2.fastq.gz" \
                        -t 'stranded_umi' \
                        -o /afs/crc/group/TIMELab/SpaceU87_Nov2024/CogentAP_OUT/ \
                        -b /afs/crc/group/TIMELab/SpaceU87_Nov2024/WellList.txt
