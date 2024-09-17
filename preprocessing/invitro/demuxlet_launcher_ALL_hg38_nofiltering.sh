#!/bin/bash

echo $LSB_JOBINDEX
echo $1
data=$(head -n $LSB_JOBINDEX $1 | tail -n1)
echo $data

INPUTBAM=$(echo $data | awk {'print $1'} )
poolnum=$(echo $data | awk {'print $2'})
INPUTVCF=genotypes/inVitroLines_hg38LO.vcf
DEMUXLET=/software/demuxlet/demuxlet
rootfolder=$(echo $INPUTBAM | sed 's/possorted_genome_bam.bam//g')
INPUTBARCODES=$(echo ${rootfolder}filtered_feature_bc_matrix/barcodes.tsv.gz)
OUTPUT=$(echo ${rootfolder}output.demuxlet.doublet0.5.noFilter)
SAMPLELIST=${poolnum}_sample_list.txt
$DEMUXLET --sam $INPUTBAM --vcf $INPUTVCF --field GT --doublet-prior 0.5 --group-list $INPUTBARCODES --out $OUTPUT --sm-list $SAMPLELIST

