#!/bin/bash
# Expect one input:
SRA=$1

# Subset alignments from Pik-1 to Pik-2
samtools view -b -h Sequencing/Alignment/N22_full/${SRA}_N22_sorted_marked.bam.cram "CM007637.2:30879574-30891650" > Variants/Vars/${SRA}_N22_sorted_marked_subset.bam

# Call variants
freebayes -f ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna Variants/Vars/${SRA}_N22_sorted_marked_subset.bam > Variants/Vars/${SRA}_N22.vcf
rm Variants/Vars/${SRA}_N22_sorted_marked_subset.bam

# Filter low quality calls
vcftools --vcf Variants/Vars/${SRA}_N22.vcf --minQ 20 --recode --recode-INFO-all --out Variants/Vars/${SRA}_N22_q20
rm Variants/Vars/${SRA}_N22.vcf

# Apply variants to reference genome
bgzip Variants/Vars/${SRA}_N22_q20.recode.vcf
tabix -p vcf Variants/Vars/${SRA}_N22_q20.recode.vcf.gz
## Pik1 Exon1
samtools faidx ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna CM007637.2:30885360-30885923 | bcftools consensus Variants/Vars/${SRA}_N22_q20.recode.vcf.gz -o Variants/Vars/temp_${SRA}_p1e1.fa
## Pik1 Exon2
samtools faidx ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna CM007637.2:30884643-30885209 | bcftools consensus Variants/Vars/${SRA}_N22_q20.recode.vcf.gz -o Variants/Vars/temp_${SRA}_p1e2.fa
## Pik1 Exon3
samtools faidx ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna CM007637.2:30879574-30881871 | bcftools consensus Variants/Vars/${SRA}_N22_q20.recode.vcf.gz -o Variants/Vars/temp_${SRA}_p1e3.fa
## Pik2 Exon1
samtools faidx ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna CM007637.2:30888421-30889410 | bcftools consensus Variants/Vars/${SRA}_N22_q20.recode.vcf.gz -o Variants/Vars/temp_${SRA}_p2e1.fa
## Pik2 Exon2
samtools faidx ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna CM007637.2:30889575-30891650 | bcftools consensus Variants/Vars/${SRA}_N22_q20.recode.vcf.gz -o Variants/Vars/temp_${SRA}_p2e2.fa

# Reverse compliment Pik1 sequences and format gene names
samtools faidx -i Variants/Vars/temp_${SRA}_p1e1.fa CM007637.2:30885360-30885923 | sed "s/\*//g" | sed "1s/.*/>${SRA}_P1E1/" > Variants/Vars/${SRA}_P1E1.fa
samtools faidx -i Variants/Vars/temp_${SRA}_p1e2.fa CM007637.2:30884643-30885209 | sed "s/\*//g" | sed "1s/.*/>${SRA}_P1E2/" > Variants/Vars/${SRA}_P1E2.fa
samtools faidx -i Variants/Vars/temp_${SRA}_p1e3.fa CM007637.2:30879574-30881871 | sed "s/\*//g" | sed "1s/.*/>${SRA}_P1E3/" > Variants/Vars/${SRA}_P1E3.fa
samtools faidx Variants/Vars/temp_${SRA}_p2e1.fa CM007637.2:30888421-30889410 | sed "s/\*//g" | sed "1s/.*/>${SRA}_P2E1/" > Variants/Vars/${SRA}_P2E1.fa
samtools faidx Variants/Vars/temp_${SRA}_p2e2.fa CM007637.2:30889575-30891650 | sed "s/\*//g" | sed "1s/.*/>${SRA}_P2E2/" > Variants/Vars/${SRA}_P2E2.fa

# Translate variants
transeq Variants/Vars/${SRA}_P1E1.fa Variants/Vars/${SRA}_P1E1.pep
transeq Variants/Vars/${SRA}_P1E2.fa Variants/Vars/${SRA}_P1E2.pep
transeq Variants/Vars/${SRA}_P1E3.fa Variants/Vars/${SRA}_P1E3.pep
transeq Variants/Vars/${SRA}_P2E1.fa Variants/Vars/${SRA}_P2E1.pep
transeq Variants/Vars/${SRA}_P2E2.fa Variants/Vars/${SRA}_P2E2.pep

# Clean up files
rm Variants/Vars/temp_${SRA}_p1e1.fa
rm Variants/Vars/temp_${SRA}_p1e1.fa.fai
rm Variants/Vars/temp_${SRA}_p1e2.fa
rm Variants/Vars/temp_${SRA}_p1e2.fa.fai
rm Variants/Vars/temp_${SRA}_p1e3.fa
rm Variants/Vars/temp_${SRA}_p1e3.fa.fai
rm Variants/Vars/temp_${SRA}_p2e1.fa
rm Variants/Vars/temp_${SRA}_p2e1.fa.fai
rm Variants/Vars/temp_${SRA}_p2e2.fa
rm Variants/Vars/temp_${SRA}_p2e2.fa.fai
