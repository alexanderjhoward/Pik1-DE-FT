# Calling *Pik-1* variants from the 3,000 Rice Genomes Project

## Background
With our fine-tuned models finished, we wanted to apply them towards the phenotyping of naturally-evolved receptor variants. Our approach was to identify *Pik-1* ligand binding domain variants within the 3k RGP dataset with complete read coverage against a reference *Pik-1* allele and use those samples for variant calling. Later on in the project we also do variant calling on the full-length *Pik-1* and *Pik-2* receptor sequences to synthesize them for agroinfiltration HR assays. 

## Method

Start by creating and activating a conda environment with all necessary dependencies.
```bash

    conda env create -f environment.yml
    conda activate 3kRGP_VC

```

We next downloaded the N22 reference rice genome which posesses a known *Pik-1* alelle (Pikp-1) that is very closely related to the allele we used in directed evolution (Pikh-1). This genome was then indexed. 

```bash

    wget -P Source/N22/  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/952/365/GCA_001952365.3_OsN22RS2/GCA_001952365.3_OsN22RS2_genomic.fna.gz
    gunzip Source/N22/GCA_001952365.3_OsN22RS2_genomic.fna.gz
    bwa index -p Sequencing/Reference_Data/N22 Source/N22/GCA_001952365.3_OsN22RS2_genomic.fna

```

Next we obtained a list of NCBI SRA IDs from [the 3k RGP BioProject](https://www.ncbi.nlm.nih.gov/biosample?Db=biosample&DbFrom=bioproject&Cmd=Link&LinkName=bioproject_biosample&LinkReadableName=BioSample&ordinalpos=1&IdsFromResult=262761) to download our sample reads. To do this, I chose "Send to" > "File" > "Accessions List" and downloaded the resulting file, named **biosample_result.txt**. After reformatting this file, I saved it as **sra.txt**, which is available in the "Source" directory already.

```bash

    dos2unix biosample_result.txt
    sed '/^$/d' biosample_result.txt > Source/sra.txt
    rm biosample_result.txt

```

To fit within computational constraints and more rapidly determine the presence/absence of Pikp-1 variants in our dataset, one lane of reads from each rice line was aligned against the N22 genome. These output alignments were checked for read coverage within the ligand binding domain, and lines with 6 or more positions missing a mapped read within that window were excluded from further consideration. Lines which passed this threshold were retained for downstream alignment using their full set of reads. This approach helped narrow the number of full sequence alignments we needed to conduct to identify variants closely related to Pikp/Pikh.

```bash

    while read i; do
        bash Scripts/bwa_pipeline_search.sh $i
    done < Source/sra.txt

```

The list of rice lines which passed the read coverage threshold was collected.

```bash

    cd Sequencing/Alignment/N22
    ls *.txt | sed 's/_N22_pileup.txt//' >> ../../../Source/N22_full_sra.txt
    cd ../../..

```

This subset of high-coverage rice varieties was re-aligned against the reference genome using the full set of reads available per line. 

```bash

    while read i; do
        bash Scripts/bwa_pipeline_full_search.sh $i
    done < Source/N22_full_sra.txt

```

Like last time, the list of successful alignments was collected. This time around, we only kept results from lines with full coverage along the ligand binding domain for variant calling.

```bash

    cd Sequencing/Alignment/N22_full
    ls *.txt | sed 's/_N22_pileup.txt//' >> ../../../Source/N22_full_vcf.txt
    cd ../../..

```

FreeBayes was used for variant calling against the reference genome and obtaining *Pik-1* sequence variants.

```bash

    while read i; do
        bash Scripts/vcf_to_fa.sh $i
    done < Source/N22_full_vcf.txt

```

After conducting our G2P anallysis, we identified two variants ("SHZ-2" and "VK") with promising enhanced Avr-PikC recognition. The SRAs of rice lines with either of these ligand binding domain variants were saved to a shortlist for variant calling, named **vars_sra.txt** which is found in the "Source" directory. This list was then used to do variant calling along the full *Pik-1* and *Pik-2* sequences using the sequence alignments we already generated.

```bash

    while read i; do
        bash Scripts/vcf_to_fa_full_seq.sh $i
    done < Source/vars_sra.txt

```

With all variant calling finished, the DNA and peptide sequences were collected into FASTA and CSV files for downstream use.

```bash

    bash Scripts/format_output.sh

```

These sequences were cleaned and assembled with the **sequence_processing.R** script. The last step was to conduct a multiple sequence alignment with MUSCLE on the **Pikp_HMA_Variants_Unique.fasta** file for downstream use.

```bash

    muscle -in Output/Pikp1_HMA_Variants_Unique.fa -out Output/Data/Pikp1_HMA_Variants_Unique_Aligned.fa
    sed 's/^>/\x00&/' Output/Pikp1_HMA_Variants_Unique_Aligned.fa | sort -z | tr -d '\0' | sed 's/^>..;/>/g' >  Output/Pikp1_HMA_Variants_Unique_Aligned_Formatted.fa
    rm Output/Pikp1_HMA_Variants_Unique_Aligned.fa

```
