# PrecisionFDA Challenge #1: Calling and Quantifying in Human RNAseq

## Introduction
I'm trying to get more bioinformatical!  As part of that, I'm hoping to work through the challenges on PrecisionFDA.  Their website claims that they aspire to be:

>A secure, collaborative, high-performance computing platform that builds a community of experts around the analysis of biological datasets in order to advance precision medicine, inform regulatory science, and enable improvements in health outcomes.

They open up challenges every 4 months or so where bioinformatical teams can compete and share their work on the topic of interest.  It seems kind of like kaggle, but with a biology bend.

You can find [their first challenge here](https://precision.fda.gov/challenges/1).  In it, they started with FASTQ sequences from Illumina reads of NA12878, a cell line from a woman from Utah.  They modified these reads to "add a number of specific SNV and InDel variants at 20% or greater variant frequency". Here's some of the information from their public website:


| <!-- -->    | <!-- -->    |
|-------------|-------------|
|Library Prep 	| Libraries were prepared from 1 microgram of genomic DNA using Kapa Biosystems' KAPA LTP Library Preparation Kit.  Libraries were enriched for exome content using Nimblegen’s SeqCap EZ Human Exome+UTR |
|Read length | 2x126bp |
|Insert size | 324bp |
|Approximate coverage | 76x |

The title of the challenge is odd as it seems to be _Whole Exome Sequencing_ rather than _RNAseq_.

Since SNP-calling is new to me, I needed a gentle introduction, and I found some helpful tips in the Biostar Handbook.  Accordingly, my first stuff here follows their examples with the Ebola genome.

## Initial plan

I'll explain these more below, but at the onset, I'm assuming that these will 
include:
	1. Practicing mapping to the Ebola virus genome (per the Biostar handbook)
		a. Mapping with real RNAseq from a paper
		b. As above, but generating fake reads and then validating that I can map them.
		c. Same fake reads as above, but ensuring I can map AND QUANTIFY them.
	2. Switching to human
		a. Mapping with subset of real reads
		b. Mapping with all reads


## Explanation of files present
The project path is probably best explained by going through the files in it.  

### Early files
Here's the earlier files run on smaller genome -- either Ebola or just the exomes of Hs chromosome 15.

- `intro_vc.sh` is a concise, simple script that can simulate some reads from a genome, intoducing some mutations into them.  It them runs a variant caller on that output.  
- `ebola_mapping.sh` is conceptually similar to the above script, except it reads in information from the command line, enabling you to choose.
	1. If you want to simulate read or download a specific set from SRA
	2. If you want to map and call variants on those reads with a `bcftools` or a `gatk` method

### Later files
- `download_fda_files.sh` downloads various files needed to do the analysis, such as the human genome, as well as the modified FASTQ files from the FDA site (note: you need permissions to access them, and the links change regularly).
- `trim_and_map.sh` is the meat of the project.  The "trimming" refers to trimming the FASTQ files to get rid of (portoins of) reads with poor quality scores.  You need to tell it through command-line argument if you're going to work with `all` of a `subset` of the reads.  It might be a to try with a subset first to make sure all the scripts run well.  I ended up using `gatk Mutect2` for predictions.  Mutect2 is for calling somatic mutations and I'll try to get into why below. 
- `run_tonight.sh` is a about 3-line script that just runs the above script, first running the subset and then running all.  I had it go overnight on AWS.
- `extract_chr15.sh` the output of `trim_and_map.sh` is LARGE and a pain to move around, so I extracted everything from just one chromosome, chr15, to download from instance and look at locally on IGV.
- `trim_mutect_outputs.sh` One thing apparent from the above was that the variant callers would happily call variants where there was only 1 read!  This was, perhaps obviously, generally reads that were far from the enriched regions.  This script filters your VCF files to only be the ones that map within 150 bp of the enriched regions.

## Selection of variant caller
I was drawn to `gatk Mutect2` because the challenge is pretty clear that it wants to know the exact rates of mutations and they won't be appearing at straightforward ratios.  The variant callers I explored in the Ebola portion were designed for calling mutations in organisms of fixed ploidy, e.g. haploid (1 of each chromosome) or diploid (2 of each chromosome).  The callers then try to shoehorn allele frequencies in to fit these expectations.   So if an allele appears 35% of the time, it might be called as 50.0% of the time, representing one of the 2 alleles in a diploid organism.  

Mutect2 avoids that by worrying about mutations in somatic cells, i.e. cells throughout the body.  It may have a lot of other weird baggage, though, so I will need to learn more about it and/or possibly try some other methods down the line, e.g. freebayes.

## Submitting to the FDA
The actual challenge ended 3 years ago, but you can still run your predictions through their app to check them.  
The first one I submitted a VCF file that had been built by Mutect2 and then trimmed by `trim_mutect_output.sh`.  It failed their quality specifications.  It was especially upset that there were too many multiallelic sites.  
I then took the VCF that had been built by Mutect2, further annotated by FilterMutectCalls, and trimmed to exome regions as above.  In order to only get higher quality calls, I ran the following on it:
```bcftools view -M 2 -f PASS filtered-bed-trimmed.vcf.gz | bgzip -c > filtered-bed-trimmed-passed.vcf.gz```
The `-M 2` limits to ones with 2 alleles only, avoiding mult-allelic issues.  That limit ultimaely doesn't matter as filtering for "PASS" already removes those.  It didn't fail this time, but it wasn't wonderful either.

|                               | Control    | Me       |
|-------------------------------|------------|----------|
| Injected indels 1\-3 found    | 11         | 9        |
| Injected indels 4\-6 found    | 11         | 11       |
| Injected indels 7\-9 found    | 11         | 9        |
| Injected indels 10\-12 found  | 11         | 7        |
| Injected indels 13\+ found    | 1          | 0        |
| Injected snps found           | 5          | 5        |
| Total injected variants found | 50         | 41       |
| Whole\-exome SNP precision    | 1          | 0\.60061 |
| Whole\-exome SNP recall       | 6\.90E\-05 | 0\.5422  |
| Whole\-exome SNP F\-score     | 0\.000139  | 0\.56991 |
| Whole\-exome indel precision  | 1          | 0\.53506 |
| Whole\-exome indel recall     | 0\.004787  | 0\.33766 |
| Whole\-exome indel F\-score   | 0\.009529  | 0\.41404 |

Obviously I missed some indels.  I'd need to dig into what the various scores on the lower half of the table mean, and think about ways to improve all of it.  Nonetheless, it was nice to get something running!

I may continue working on this, but a new challenge (#12) is opening in a day or two (2020-11-30), and I think I'll want to move onto that and see if I can get good enough at that to throw my hat in the ring for a real running one!