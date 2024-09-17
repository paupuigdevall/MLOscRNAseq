#!/bin/bash

echo $LSB_JOBINDEX
echo $1
line=$(head -n $LSB_JOBINDEX $1 | tail -n1)
echo $line

sampleid=$(echo $line | awk {'print $1'})
folder=$(echo $line | awk {'print $2'})

TPATH=refgenomes/cellranger/refdata-cellranger-GRCh38-3.0.0

/software/kilpinen/cellranger-3.1.0/cellranger count --id=${sampleid} --fastqs=${folder} --transcriptome=${TPATH} --jobmode=local --localcores=32 --localmem=200 --maxjobs=64
