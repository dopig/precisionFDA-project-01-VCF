#!/bin/bash


#
# It's expensive to pull all of the many files from AWS to my laptop in order
# to inspect them on IGV.  Also, my computer doesn't like keeping that much 
# in memory.  Therefore, I'm going to just pull parts of all the files from
# one chromosome, number 15, into a zip file that I'll download and inspect
# locally.
#

set -euo pipefail


INDIR=results/all

BAM1=align_initial.bam
BAM2=align_marked_dups.bam
BAM3=align_marked_dups_bqsr.bam
VCF1=unfiltered.vcf
VCF2=filtered.vcf

OUTDIR=${INDIR}/all-chr15


# First make a chromosome-15 directory

if ! test -d "$OUTDIR"; then
    echo "Making directory $OUTDIR"
    mkdir -p $OUTDIR
fi


# Then deal with the BAM files, by pulling out chromosome 15 & by reindexing them

for file in $BAM1 $BAM2 $BAM3; do
	samtools view ${INDIR}/${file} 15 -hb > ${OUTDIR}/chr15-${file}
	samtools index ${OUTDIR}/chr15-${file}
done


# Then deal with the VCF files, by pulling out chromosome 15 & by reindexing them

for file in $VCF1 $VCF2; do
	bcftools view ${INDIR}/${file} --targets 15 -Oz -o ${OUTDIR}/chr15-${file}.gz
	bcftools index ${OUTDIR}/chr15-${file}.gz
done

zip -r ${OUTDIR}.zip $OUTDIR