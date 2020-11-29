#!/bin/bash

set -euo pipefail
set -x #echo on


ORIGINAL_DIR="refs/original"
WORKING_DIR="refs/working"
RUNLOG="refs/LOG.md"

for dir in $ORIGINAL_DIR $WORKING_DIR; do
	if ! [ -d $dir ]; then mkdir -p $dir; fi
done

# echo "*** Preparing References ***"  >> $RUNLOG
# echo "`date`: Downloading human_g1k_v37.fasta, NA12878MOD_1.fastq, & NA12878MOD_1.fastq to ${ORIGINAL_DIR}" >> $RUNLOG

# wget --directory-prefix=$ORIGINAL_DIR ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
# wget --directory-prefix=$ORIGINAL_DIR https://dl.dnanex.us/F/D/p6Bg2YGVVjbgx28kp77153G3bv3265p0BjbjF0fx/NA12878MOD_1.fastq.gz
# wget --directory-prefix=$ORIGINAL_DIR https://dl.dnanex.us/F/D/0250qP2GQyQjyYB1GkgXxPbYyGkqZ4BX6px5X327/NA12878MOD_2.fastq.gz

# echo "`date`: Making bgzip versions of those 3 files in $WORKING_DIR" >> $RUNLOG

set +euo pipefail
for filename in ${ORIGINAL_DIR}/human_g1k_v37.fasta.gz; do
	zcat ${filename} | bgzip -c \
		> ${WORKING_DIR}/$(basename $filename) \
		2>> $RUNLOG
done
set -euo pipefail

echo "`date`: Download and unzip dbsnp132_20101103.vcf to be used with 'gatk BaseRecalibrator' later" >> $RUNLOG
wget --directory-prefix=$ORIGINAL_DIR ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/dbsnp132_20101103.vcf.gz
zcat ${ORIGINAL_DIR}/dbsnp132_20101103.vcf.gz > ${WORKING_DIR}/dbsnp132_20101103.vcf  2>> $RUNLOG
echo "*** Reference Prep Done ***"  >> $RUNLOG