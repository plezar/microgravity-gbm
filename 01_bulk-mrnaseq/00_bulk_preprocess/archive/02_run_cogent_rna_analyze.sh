#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m abe            # Send mail when job begins, ends and aborts
#$ -pe smp 24     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N cogentAP       # Specify job name

$COGENT_AP_HOME/cogent rna analyze -g 'hg38' \
			-o /afs/crc/group/TIMELab/SpaceU87_Nov2024/CogentAP_OUT \
			-i /afs/crc/group/TIMELab/SpaceU87_Nov2024/ILMN_2346_HarkerBrent_ND_1Pool_Nov2024 \
			-t 'stranded_umi' \
			--ribodepletion 'auto'
