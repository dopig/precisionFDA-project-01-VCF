#!/bin/bash

set -euo pipefail


# The reference file, which is a bgzipped version of file downloaded to
# refs/original, decompressed, and recompressed with bgzip
# by download_fda_files.sh 
REF=refs/working/human_g1k_v37.fasta.gz

# The FASTQ read files are downloaded by the same script as the reference.
R1=refs/original/NA12878MOD_1.fastq.gz
R2=refs/original/NA12878MOD_2.fastq.gz

# If decide to test just a subset of the reads, this is where they'll be
R1_lil=refs/working/NA12878MOD_10000_1.fastq.gz
R2_lil=refs/working/NA12878MOD_10000_2.fastq.gz

# One of the steps leading up to gatk Mutect2 depends on this file.
# It doesn't like it to be zipped.
KNOWNSITES=refs/working/dbsnp132_20101103.vcf




# These are ones you may want to alter depending on run!

if [ $1 = "subset" ]; then
	USE_LIL=true
	OUTDIR="results/subset-of-reads"
elif [ $1 = "all" ]; then
	USE_LIL=false
	OUTDIR="results/all-reads"
else
	echo "Need to specify if operating on subset of reads or all reads."
	exit
fi


# Get some GATK stuff sorted out first


if ! test -f "${REF/.fasta*/.dict}"; then
	gatk CreateSequenceDictionary -REFERENCE $REF -OUTPUT "${REF/.fasta*/.dict}"
fi

if ! [ -f ${KNOWNSITES}.idx ]; then
	gatk IndexFeatureFile --input $KNOWNSITES
fi

if [ $USE_LIL = true ] ; then
	if ! [ -f $R2_lil ] ; then
	    seqtk sample -s85 $R1 10000 | gzip -c > $R1_lil
	    seqtk sample -s85 $R2 10000 | gzip -c > $R2_lil
	fi
	R1=$R1_lil
	R2=$R2_lil
fi


# The name of the BAM file
BAM=${OUTDIR}/align_initial.bam
BAM_DUP=${OUTDIR}/align_marked_dups.bam
BAM_BQSR=${OUTDIR}/align_marked_dups_bqsr.bam

# Make outdir
if [  true = false ] ; then
	if ! test -d "$OUTDIR"; then
	    echo "Making directory $OUTDIR"
	    mkdir -p $OUTDIR
	else
		echo "This out-directory already exists."
		exit
	fi
fi


# Create a bwa index for the reference.
if ! test -f "${REF}.sa"; then
    echo "Creating bwa index for $REF"
    bwa index $REF
fi


# This file will keep track of messages printed during the run.
RUNLOG=${OUTDIR}/runlog.txt
echo "*** Run by `whoami` on `date`" > $RUNLOG
echo "*** Run with $R1 and $R2 against $REF" >> $RUNLOG




if [  true = false ] ; then

	echo "Trimming with bbduk.sh at `date`" >> $RUNLOG

	BBDUK1=${OUTDIR}/passed_bbduk_1.fq
	BBDUK2=${OUTDIR}/passed_bbduk_2.fq
	BBDUK1_failed=${OUTDIR}/failed_bbduk_1.fq
	BBDUK2_failed=${OUTDIR}/failed_bbduk_2.fq


	bbduk.sh in1=$R1 in2=$R2 ref=adapters \
		outm1=$BBDUK1_failed out1=${BBDUK1} \
		outm2=$BBDUK2_failed out2=${BBDUK2} \
		qtrim=r overwrite=true qtrim=30 mlf=0.5 2>> $RUNLOG



	# I didn't have this tag when I did it the first time, and it took hours. :-/

	if ! test -f "${BAM}.bai"; then
		echo "bwa aligning at `date`" >> $RUNLOG
		TAG="@RG\tID:precisionFDA\tSM:precisionFDA-1\tLB:precisionFDA\tPL:ILLUMINA"	# GATK requires a read-group header in $BAM file
		bwa mem -R $TAG $REF $BBDUK1 $BBDUK2 | samtools sort > $BAM 2>> $RUNLOG
		#bwa mem $REF $R1 $R2 | samtools sort > $BAM 2>> $RUNLOG

		# Index the BAM file
		echo "indexing ${BAM} at `date`" >> $RUNLOG
		samtools index $BAM 2>> $RUNLOG
	fi


	echo ">>>>> Now that the BAM alignment has been made, note and delete the bbduk files <<<<<" >> $RUNLOG

	for filename in $BBDUK1 $BBDUK2 $BBDUK1_failed $BBDUK2_failed; do
		echo `ls -lrth $filename`
		rm $filename
	done
fi

# As here: https://gatk.broadinstitute.org/hc/en-us/articles/360050814112-MarkDuplicatesSpark
echo "*** Time: `date` : Running gatk MarkDuplicatesSpark" >> $RUNLOG
gatk MarkDuplicatesSpark \
    -I $BAM \
    -O $BAM_DUP \
    -M ${OUTDIR}/marked_dup_metrics.txt 2>> $RUNLOG
    #--remove-sequencing-duplicates


echo "*** Time: `date` : Running gatk BaseRecalibrator" >> $RUNLOG
gatk BaseRecalibrator \
	-I $BAM_DUP \
	-R $REF \
	--known-sites $KNOWNSITES \
	-O ${OUTDIR}/recal_data.table 2>> $RUNLOG

echo "*** Time: `date` : Running gatk ApplyBQSR" >> $RUNLOG
gatk ApplyBQSR \
	-R $REF \
	-I $BAM_DUP \
	--bqsr-recal-file ${OUTDIR}/recal_data.table \
	-O $BAM_BQSR 2>> $RUNLOG


gatk Mutect2 -R $REF -I $BAM_BQSR -O ${OUTDIR}/unfiltered.vcf

# Optional-1: GetPileupSummaries, CalculateContamination
# Optional-2: LearnReadOrientationModel --> this is especially if FFPE samples

gatk FilterMutectCalls -R $REF -V ${OUTDIR}/unfiltered.vcf -O ${OUTDIR}/filtered.vcf
