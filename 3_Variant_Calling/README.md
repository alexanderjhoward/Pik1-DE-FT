# Calling *Pik-1* variants from the 3,000 Rice Genomes Project

## Background
With our fine-tuned models finished, we wanted to apply them towards the phenotyping of naturally-evolved receptor variants. Out approach was to identify *Pik-1* ligand binding domain variants within the 3,000 Rice Genomes Project dataset with complete read coverage in the *Pik-1* ligand binding domain and use those samples for variant calling. Later on in the project we also do variant calling on the full-length *Pik-1* and *Pik-2* receptors to synthesize them for agroinfiltration HR assays. 

## Method
I used BWA for read alignment. I started by indexing our reference genome N22, which posesses a known *Pik-1* alelle (*Pikp-1/2*) that is closely related to the allele we used in directed evolution (*Pikh-1*).

 ```bash

    bwa index -p Sequencing/Reference_Data/N22 ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna

```

Next we obtained a list of NCBI IDs to download our sample reads from. I went to the 3k RGP BioProject and [retreived all listed SRAs](https://www.ncbi.nlm.nih.gov/biosample?Db=biosample&DbFrom=bioproject&Cmd=Link&LinkName=bioproject_biosample&LinkReadableName=BioSample&ordinalpos=1&IdsFromResult=262761). To do this, I chose "Send to" > "File" > "Accessions List" and downloaded the resulting file, named **biosample_result.txt**. After reformatting this file, I saved it as **sra.txt**, which is available in the "Source" directory already.

```bash

    dos2unix biosample_result.txt
    sed '/^$/d' biosample_result.txt > Source/sra.txt
    rm biosample_result.txt

```

To rapidly determine the presence/absence of *Pikp-1* variants in our dataset, one lane of reads from each rice line was aligned against the N22 genome. Output alignments were checked for read coverage within the ligand binding domain, and lines with 6 or more positions missing a read within that window were excluded from further consideration. Lines with full read coverage or missing 5 or fewer positions were retained for alignment with that line's full set of reads.

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
        bash Scripts/bwa_pipeline_n22_full.sh $i
    done < Source/N22_full_sra.txt

```

Like last time, the list of successful alignments was collected. This time around, we only kept results from lines with full coverage along the ligand binding domain for variant calling.

```bash

    cd Sequencing/Alignment/N22_full
    ls *.txt | sed 's/_N22_pileup.txt//' >> ../../../Source/N22_full_vcf.txt
    cd ../../..

```

FreeBayes was used to do variant calling against the reference genome and obtain *Pik-1* sequence variants.

```bash

    while read i; do
        bash Scripts/N22_vcf_to_fa.sh $i
    done < Source/N22_full_vcf.txt

```

After conducting our G2P anallysis, two variants ("SHZ-2" and "VK") which showed promise for enhanced ligand recognition. The SRAs of rice lines with either of these ligand binding domain variants were saved to a shortlist for variant calling, named **vars_sra.txt** which is found in the "Source" directory. This list was then used to do variant calling on the preexisting sequence alignments we generated already, just along the full *Pik-1* and *Pik-2* sequences.

```bash

    while read i; do
        bash Scripts/Vars_vcf_to_fa.sh $i
    done < Source/vars_sra.txt

```

With all variant calling finished, the DNA and peptide sequences were collected into FASTA and CSV files for downstream use.

```bash

    bash Scripts/format_output.sh

```

These sequences were cleaned and assembled with the **sequence_processing.R** script. The last step taken was to conduct a multiple sequence alignment with MUSCLE on the **Pikp_HMA_Variants_Unique.fasta** file for downstream use.

```bash

    muscle -in Output/Pikp1_HMA_Variants_Unique.fa -out Output/Data/Pikp1_HMA_Variants_Unique_Aligned.fa
    sed 's/^>/\x00&/' Output/Pikp1_HMA_Variants_Unique_Aligned.fa | sort -z | tr -d '\0' | sed 's/^>..;/>/g' >  Output/Pikp1_HMA_Variants_Unique_Aligned_Formatted.fa
    rm Output/Pikp1_HMA_Variants_Unique_Aligned.fa

```
