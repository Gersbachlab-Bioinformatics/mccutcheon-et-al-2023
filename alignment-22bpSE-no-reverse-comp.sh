#! /bin/bash

## command line parameters
# $1: Reference genome location
# $2: Folder with FASTQ files
# $3: Sample Read 1 FASTQ file
# $4: Output basename
# $5: Protospacers FASTA file [Optional, only required the first time]

## dependencies: 
# - samtools >=1.3.1
# - bowtie2 >=2.3.5.1

# Build Bowtie2 custom reference index if it doesn't exist
if [[ ! -e $1.1.bt2 ]];
then
	bowtie2-build $5 $1
fi


perform_alignment ()
{

# Map 22bp reads in end-to-end --very-sensitive mode 
bowtie2 \
	--very-sensitive \
	--end-to-end \
	--threads 32 `# Using 32 CPUs, adjust as needed` \
	--norc `# Do not try to map reverse complements` \
	-3 1 `#trimming 1bp from 3' (as recommended by sequencing core)` \
	-I 0 `# no mimimim length alignment` \
	-X 200 `# to reveal chimeric alignments, but faster than the 500bp default` \
	-x $1 `# Bowtie2 index basename` \
	-U $2/$3 `# FASTQ file with reads to align` \
| samtools view -bS - > $4.bam `# Save BAM file` \
&& samtools sort \
	-@ 32 `# Using 32 CPUs, adjust as needed` \
	-m 100G  `# Using 100G of memory, adjust as needed` \
	$4.bam \
	-o $4.sorted.bam  \
&& samtools index $4.sorted.bam \
&& samtools idxstats `# Extract count table for each element (protospacer) in the reference` \
	$4.sorted.bam \
	| awk '{print $1"\t"$3}' \
	| grep -v "^\*" \
> $4.counts.txt

# Create shuffled BAM file for library complexity estimation
samtools bamshuf $4.sorted.bam $4.shuffle 

}

perform_alignment $1 $2 $3 $4
