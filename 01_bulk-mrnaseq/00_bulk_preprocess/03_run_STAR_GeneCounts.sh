#!/bin/bash

#$ -M mzarodn2@nd.edu   # Email address for job notification
#$ -m ae            # Send mail when job begins, ends and aborts
#$ -pe smp 8     # Specify parallel environment and legal core size
#$ -q long           # Specify queue
#$ -N STAR       # Specify job name

# this remains the same for all jobs (specific to mapping)
OUTPREFIX=$1
NCPU=12
INPUTDIR=$2
OUTDIR=$3
GTFFILE=$4
GENOMEDIR=$5
INPUTDIR_trimmed="/afs/crc/group/TIMELab/SpaceU87_Nov2024/UMITOOLS_out/trimmed/"

module load bio/star/2.7.2

mkdir -p $OUTDIR/$OUTPREFIX
cd $OUTDIR/$OUTPREFIX

echo "Sample: ${OUTPREFIX}"
echo "R1 path: ${INPUTDIR}/${OUTPREFIX}_R1_001.fastq.gz"
echo "R2 path: ${INPUTDIR_trimmed}/${OUTPREFIX}_R2_001.fastq.gz"

STAR \
	--outSAMattributes All \
	--outSAMtype BAM SortedByCoordinate \
	--runThreadN $NCPU \
	--sjdbGTFfile $GTFFILE \
	--outMultimapperOrder Random \
	--genomeDir $GENOMEDIR \
	--readFilesIn <(gunzip -c ${INPUTDIR}/${OUTPREFIX}_R1_001.fastq.gz) <(gunzip -c ${INPUTDIR_trimmed}/${OUTPREFIX}_R2_001.fastq.gz) \
	--outFileNamePrefix $OUTPREFIX

