#!/bin/bash -login
# Expect one input:
SRA=$1

# Subset alignments
samtools view -b -h Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam.cram "CM007637.2:30884985-30885368" > Variants/N22/${SRA}_N22_sorted_marked_subset.bam

# Call variants
freebayes -f ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna Variants/N22/${SRA}_N22_sorted_marked_subset.bam > Variants/N22/${SRA}_N22.vcf
rm Variants/N22/${SRA}_N22_sorted_marked_subset.bam

# Filter low quality calls
vcftools --vcf Variants/N22/${SRA}_N22.vcf --minQ 20 --recode --recode-INFO-all --out Variants/N22/${SRA}_N22_q20
rm Variants/N22/${SRA}_N22.vcf

# Apply variants to reference genome
bgzip Variants/N22/${SRA}_N22_q20.recode.vcf
tabix -p vcf Variants/N22/${SRA}_N22_q20.recode.vcf.gz
samtools faidx ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna CM007637.2:30885360-30885368 | bcftools consensus Variants/N22/${SRA}_N22_q20.recode.vcf.gz -o Variants/N22/temp_${SRA}_e1.fa
samtools faidx ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna CM007637.2:30884985-30885209 | bcftools consensus Variants/N22/${SRA}_N22_q20.recode.vcf.gz -o Variants/N22/temp_${SRA}_e2.fa

# Reverse compliment variant sequences and format gene names
samtools faidx -i Variants/N22/temp_${SRA}_e1.fa CM007637.2:30885360-30885368 | sed "s/\*//g" | sed "1s/.*/>${SRA}_E1/" > Variants/N22/${SRA}_E1.fa
samtools faidx -i Variants/N22/temp_${SRA}_e2.fa CM007637.2:30884985-30885209 | sed "s/\*//g" | sed "1s/.*/>${SRA}_E2/" > Variants/N22/${SRA}_E2.fa

# Translate variants
transeq Variants/N22/${SRA}_E1.fa Variants/N22/${SRA}_E1.pep
transeq Variants/N22/${SRA}_E2.fa Variants/N22/${SRA}_E2.pep

# Clean up files
rm Variants/N22/temp_${SRA}_e1.fa
rm Variants/N22/temp_${SRA}_e1.fa.fai
rm Variants/N22/temp_${SRA}_e2.fa
rm Variants/N22/temp_${SRA}_e2.fa.fai
