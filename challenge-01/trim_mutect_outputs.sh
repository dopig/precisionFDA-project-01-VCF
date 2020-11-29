#!/bin/bash

set -euo pipefail


# The goal of this is to trim existing GATK outputs.
# Inputs: Existing Mutect2 and FilterMutectCall VCF files
# They will be trimmed based on being within 150 bp of the exome enrichments.
#
# Of course this trimming is something that can be done during the initial
# variant calling, if I were to redo it.  
# 
# I was initially going to rerun everything using the Panel of Normals (PON) et 
# al, but it seems more complicated than it should be to pull things out of 
# Google Cloud Buckets from command line, so I will just trim what I already 
# have.

#
# Get relevant kit's bed file laying out where the enriched genome regions are.
#
ORIGINAL_BED=refs/original/SeqCap_EZ_ExomeV3_Plus_UTR_hg19_UCSCBrowser.bed

if ! [ -f $ORIGINAL_BED ]; then
	curl https://sequencing.roche.com/content/dam/rochesequence/worldwide/shared-designs/Exome_UTR_Design_Annotation_Files.zip > temporary.zip
	unzip -d refs/original temporary.zip $(basename $ORIGINAL_BED)
	rm temporary.zip
fi


#
# Pad the BED file by 150 nt on each side as it seems to have a lot of coverage
# adjacent to select regions.
#
# Bedtools says it needs to to know how big each chromosome is so it doesn't
# pad too far.  I'd be content with overpadding here, but I don't think it 
# will run without this information, so I'll generate that information here
# as $GENOME_TXT.
# 

REF=refs/working/human_g1k_v37.fasta.gz
GENOME_TXT=refs/working/human_g1k_v37.txt  
# Gather the name and length, removing scaffolds after the MT (which is #25)
cut -f 1,2 ${REF}.fai | head -25 > $GENOME_TXT

PADDED_BED=refs/working/SeqCap_EZ_ExomeV3_Plus_UTR_hg19_UCSCBrowser_padded_150.bed

bedtools slop \
	-i $ORIGINAL_BED \
	-g $GENOME_TXT \
	-b 150 \
	| sort -k1,1 -k2,2n \
	| bedtools merge \
	| sed 's/chr//' \
	> $PADDED_BED


for root in results/trimmed_all_a/unfiltered results/trimmed_all_a/filtered; do

	if ! [ -f ${root}.vcf.gz.tbi ]; then
		bgzip ${root}.vcf
		tabix -p vcf ${root}.vcf.gz
	fi

	bcftools view \
		--targets-file $PADDED_BED \
		-Oz -o ${root}-bed-trimmed.vcf.gz \
		${root}.vcf.gz 
done


# # Getting these files proved to be too annoying
# curl https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf | bgzip mep.vcf.gz

# "https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf"
# "https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.idx"

# # Make these cones compressed if can
# "https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf"
# "https://console.cloud.google.com/storage/browser/_details/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx"

