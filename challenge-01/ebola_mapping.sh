#!/bin/bash

# 
#  Reading "The Biostar Handbook", they're into testing Variant Calling 
#  on the nice, little Ebola genome.  They download some reads from SRA
#  and others are generated with simulators.  
#
#  This script uses some of the tools they mention (especially in chapter 88) 
#  to work with existing SRA reads.
#
#  It takes those reads, 'bwa mem' aligns them and then tries to call variants.
# 
#  Mainly trying to get things installed and make sure they run ok, and can
#  be veiwed nicely in IGV.  Not trying to make _the best calls of all time_
#  here.
#


set -euo pipefail


ACC=AF086833		# Reference accession numbers.
REF_FA=refs/ebola/${ACC}.fa
REF=${REF_FA}.gz  # The name of the reference.

mkdir -p $(dirname $REF)

#
# Ensure proper arguments supplied ('sra' vs 'sim')
#

if ([ "$#" = "0" ]) || ! ([ $1 = 'sra' ] || [ $1 = 'sim' ]); then
    explanation='In order for the script to know if you want to use *SRA*
    reads or to (make and) use *simulated* reads, you must call it as follows:'
    echo -e $explanation "\n\t$0 sra\n\t$0 sim"
    exit 1
fi


#
# Download and index reference (if samtools-indexed ref file doesn't exist)
#

if ! test -f "${REF}.fai"; then
	echo "Downloading $REF"
	efetch -db=nuccore -format=fasta -id=$ACC > $REF_FA
	bgzip $REF_FA
    echo "Creating bwa index for $REF"
    bwa index $REF 		# This is for aligner
    samtools faidx $REF # This is for IGV; (btw )faidx is "FastA InDeX")
fi


#
# Download or generate reads
#


if [ $1 = "sra" ] ; then
	echo "*Working with reads from SRA*"
	
	OUTDIR=results/srr-align
	ID=SRR1553500
	PREFIX=${OUTDIR}/${ID}
	R1=${PREFIX}_1.fastq
	R2=${PREFIX}_2.fastq
	BAM=${PREFIX}-aligned.bam
	PLOIDY=1

	#
	# Download SRR files (only if not already done)
	#

	if ! test -f "${R2}.gz"; then 
		mkdir -p $OUTDIR
		fastq-dump -X 100000 --split-files -O $OUTDIR $ID
		for i in $R1 $R2; do bgzip $i; done
	fi

elif [ $1 = "sim" ] ; then
	echo "*Working with simulated reads*"

	OUTDIR=results/simulated-align
	ID=simulated
	PREFIX=${OUTDIR}/${ID}
	R1=${PREFIX}.bwa.read1.fastq
	R2=${PREFIX}.bwa.read2.fastq
	BAM=${PREFIX}-aligned.bam
	PLOIDY=2 # If change this to PLOIDY=1 need to add "-H" arg to dwgsim call 

	#
	# Generate simulated reads
	#

	if ! test -f "${R2}.gz"; then 
		message="Unfortunately dwgsim will not operate on compressed genomes,
			so need to unzip, simulate, and then delete."
		gunzip < $REF > $REF_FA

		# Simulate reads from the reference file.
		# dwgsim -H -r 0 -R 0.5 -X 0.7 $REF ${OUTDIR}/simulated  
		# For NO SNPs -> `-r 0`
		mkdir -p $OUTDIR
		dwgsim \
			-F 0.8 \
			-R 0.5 \
			-X 0.7 \
			-1 126 -2 126 -d 576 \
			-C 76 \
			$REF_FA $PREFIX
		rm $REF_FA
		rm ${PREFIX}.bfast.fastq
		for i in $R1 $R2; do bgzip $i; done
	fi

fi

#
# Generate the alignment (if not already done)
#

if ! test -f "${BAM}.bai"; then
	echo "*Aligning reads with bwa mem*"
	# GATK requires a read-group header in $BAM file
	TAG="@RG\tID:${ID}\tSM:${ID}\tLB:${ID}"	
	bwa mem -R $TAG $REF $R1.gz $R2.gz | samtools sort > $BAM

	# Index the BAM file
	samtools index $BAM
fi


#
# Call variants (if `bcftools` or `gatk` supplied as command-line argument)
#

if [ "$#" -le 1 ] || ! ([ $2 = 'bcftools' ] || [ $2 = 'gatk' ]); then
    explanation='If you would like to proceed to variant calling with the 
	alignments you just made, add either *bcftools* or *gatk* to your arguments 
	as follows:'
    echo -e $explanation "\n\t$0 $1 bcftools\n\t$0 $1 gatk"
    exit 0

elif [ $2 = "bcftools" ]; then
	echo "running bcftools"

	# Frustratingly, you can tell bcftools call "-ploidy 1" for haploid calls,
	# but if you tel it "-ploidy 2", which is it's default (!!) it gets angry,
	# so here's my workaround.

	if [ $PLOIDY -eq 1 ]; then
		PLOIDY_FLAG='--ploidy 1'
	elif [ $PLOIDY -eq 2 ]; then
		PLOIDY_FLAG=''
	else 
		echo "Ploidy of '$PLOIDY' is messing things up.  Aborting."
		exit 1
	fi

	# Combine the above into one step!

	bcftools mpileup -Ovu -f $REF $BAM \
		| bcftools call $PLOIDY_FLAG -vm -Oz \
		> $OUTDIR/variants_bcf.vcf.gz

elif [ $2 = "gatk" ]; then
	echo "running gatk HaplotypeCaller"
	# Already ran `samtools faidx $REF` above
	set -x

	REF_DICT=${REF/.fa.gz/.dict}

	if ! [ -f $REF_DICT ]; then
		gatk CreateSequenceDictionary -REFERENCE $REF -OUTPUT $REF_DICT
	fi

	# For SRA, GATK takes 13 minutes on my laptop vs <13 seconds for bcftools :/
	gatk HaplotypeCaller \
		-R $REF \
		-I $BAM \
		-ploidy $PLOIDY \
		-O $OUTDIR/variants_gatk.vcf.gz
fi





