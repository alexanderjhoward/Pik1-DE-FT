#!/bin/bash
# Input sample
SRA=$1

# Save prefetch data to a temporary directory
mkdir ${SRA}_temp
prefetch --output-directory ${SRA}_temp/ ${SRA}

# Loop over prefetch data with fasterq-dump
for dir in ${SRA}_temp/*/; do
        NAME=$(echo $dir | awk -F '/' '{print $2}')
        fasterq-dump ${SRA}_temp/${NAME}/${NAME}.sra -e 8 -O ${SRA}_temp/
done

# Combine all FWD and REV reads together
cat ${SRA}_temp/*_1.fastq > ${SRA}_1.fastq
cat ${SRA}_temp/*_2.fastq > ${SRA}_2.fastq
rm -r ${SRA}_temp

# Align reads against reference with BWA-MEM
bwa mem -M -t 8 Sequencing/Reference_Data/N22 ${SRA}_1.fastq ${SRA}_2.fastq 2> Logs/BWA_FULL/${SRA}_N22_full_bwa.err > Sequencing/Alignment/N22_full/${SRA}_N22.sam
rm ${SRA}_1.fastq
rm ${SRA}_2.fastq

# Sort SAM file by coordinates
picard SortSam -I Sequencing/Alignment/N22_full/${SRA}_N22.sam -O Sequencing/Alignment/N22_full/${SRA}_N22_sorted.sam -SO coordinate --VALIDATION_STRINGENCY SILENT
rm Sequencing/Alignment/N22_full/${SRA}_N22.sam

# Mark duplicate reads and generate BAM file
picard MarkDuplicates -I Sequencing/Alignment/N22_full/${SRA}_N22_sorted.sam -O Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam -M Logs/BWA_FULL/${SRA}_N22_full_metrics.txt -AS true --VALIDATION_STRINGENCY SILENT
rm Sequencing/Alignment/N22_full/${SRA}_N22_sorted.sam

# Convert BAM to CRAM to save
samtools view -T Source/GCA_001952365.3_OsN22RS2_genomic.fna -C -o Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam.cram Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam
rm Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam

# Index CRAM file
samtools index Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam.cram

# Explore read coverage
samtools mpileup -f Source/GCA_001952365.3_OsN22RS2_genomic.fna -r CM007637.2:30884985-30885209 -s -aa Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam.cram > Sequencing/Alignment/N22_full/${SRA}_N22_pileup.txt
samtools mpileup -f Source/GCA_001952365.3_OsN22RS2_genomic.fna -r CM007637.2:30885360-30885368 -s -aa Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam.cram >> Sequencing/Alignment/N22_full/${SRA}_N22_pileup.txt

# Sort out low coverage assemblies (regions of pileup missing read alignment)
# If a given sample has any unaligned bases, move results to a separate folder
N22_Gap=$(cat Sequencing/Alignment/N22_full/${SRA}_N22_pileup.txt | awk '{ if ($4 < 1) { print } }' | wc -l)
if [ $N22_Gap -gt 0 ]
then
	mv Sequencing/Alignment/N22_full/${SRA}_N22_pileup.txt Sequencing/Alignment/N22_full/Failed/${SRA}_N22_pileup.txt
	mv Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam.cram Sequencing/Alignment/N22_full/Failed/${SRA}_N22_sorted_marked.bam.cram
	mv Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam.cram.crai Sequencing/Alignment/N22_full/Failed/${SRA}_N22_sorted_marked.bam.cram.crai
fi
