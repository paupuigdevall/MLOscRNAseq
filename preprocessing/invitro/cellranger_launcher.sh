#!/bin/bash

echo $LSB_JOBINDEX
echo $1
folder=$(head -n $LSB_JOBINDEX $1 | tail -n1)
echo $folder

TPATH=refgenomes/cellranger/refdata-cellranger-GRCh38-3.0.0
sampleid=cellranger-hg38

/software/cellranger-3.1.0/cellranger count --id=${sampleid} --fastqs=${folder} --transcriptome=${TPATH} --jobmode=local --localcores=32 --localmem=200 --maxjobs=64

