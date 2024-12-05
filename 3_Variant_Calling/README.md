# Aligning and calling *Pik-1* variants from the 3,000 Rice Genomes Project

## Background
While variants have been called on the 3k RGP before, the alignments were either not easily accessible or the reference genome lacked an intact *Pik-1* gene. The following approach was taken to identify variants within the dataset which posessed a *Pik-1* gene by doing an initial alignment with a subset of sample reads before doing a full alignment with the rice lines which showed complete or nearly complete read coverage at the *Pik-1* ligand binding domain. 

## Set up
I used BWA for read alignment. I started by indexing our reference genomes N22 (*Pikp-1/2*) and Kitaake (*Pikm-1/2*).

 ```bash

    bwa index -p Sequencing/Reference_Data/N22 ../1_Create_Databases/GenomeDB/N22/GCA_001952365.3_OsN22RS2_genomic.fna
    bwa index -p Sequencing/Reference_Data/Kit ../1_Create_Databases/GenomeDB/Kitaake/GCA_009797565.1_Osativa_Kitaake_v2.0_genomic.fna

```

Next we needed a list of NCBI IDs to download our sample reads from. I got mine from the 3k RGP BioProject and [retreiving all listed SRAs](https://www.ncbi.nlm.nih.gov/biosample?Db=biosample&DbFrom=bioproject&Cmd=Link&LinkName=bioproject_biosample&LinkReadableName=BioSample&ordinalpos=1&IdsFromResult=262761). To do this, I chose "Send to" > "File" > "Accessions List" and downloaded the resulting file, named **biosample_result.txt**. This had to be formatted into a list without blank rows (and put into Unix format in my case, since I work on Windows).

```bash

    dos2unix biosample_result.txt
    sed '/^$/d' biosample_result.txt > Source/sra.txt
    rm biosample_result.txt

```

## Initial *Pik-1* alignment
To quickly determine the presence/absence of rice varieties that posessed variants of *Pikp-1* or *Pikm-1*, a subset of reads from each rice line was aligned against both reference genomes. Output alignments were checked for read coverage within the HMA ligand binding domain, and lines that had 6 or more positions missing a read within that stretch were removed from further consideration. Any lines with full read coverage or missing 5 or fewer positions were retained for a more full-coverage alignment with all of that line's reads.

```bash

    while read i; do
        bash bwa_pipeline_search.sh $i
    done < Source/sra.txt

```

The list of rice lines with acceptable read coverage were extracted for use in next step's full alignment.

```bash

    cd Sequencing/Alignment/N22
    ls *.txt | sed 's/_N22_pileup.txt//' >> ../../../Source/N22_full_sra.txt
    cd ../Kit
    ls *.txt | sed 's/_Kit_pileup.txt//' >> ../../../Source/Kit_full_sra.txt
    cd ../..

```

## Full *Pik-1* alignment
With our subset of well-aligned rice varieties we can re-align against the relevent reference genomes using the full set of reads available per line. 

```bash

    while read i; do
        bash bwa_pipeline_n22_full.sh $i
    done < Source/N22_full_sra.txt

    while read i; do
        bash bwa_pipeline_kit_full.sh $i
    done < Source/Kit_full_sra.txt

```

Like last time, we need to collect the successful alignments. This time around, the script only keeps results from lines with full coverage along the HMA domain, so we'll just collect those for variant calling.

```bash

    cd Sequencing/Alignment/N22_full
    ls *.txt | sed 's/_N22_pileup.txt//' >> ../../../Source/N22_full_vcf.txt
    cd ../Kit_full
    ls *.txt | sed 's/_Kit_pileup.txt//' >> ../../../Source/Kit_full_vcf.txt
    cd ../..

```

## *Pik-1* ligand binding domain variant calling
After full-scale alignment of our samples of interest, we use FreeBayes to do variant calling against the reference genome to obtain the sequence variants found in the dataset.

```bash

    while read i; do
        bash N22_vcf_to_fa.sh $i
    done < Source/N22_full_vcf.txt

    while read i; do
        bash Kit_vcf_to_fa.sh $i
    done < Source/Kit_full_vcf.txt

```

## Full-length *Pik-1* variant calling
After doing machine learning, we identified two HMA variants ("San" and "Vel") which showed the most promise at enhanced ligand recognition and wanted to get the full sequences for. The SRAs of rice lines with either of these HMA domains were saved to a shortlist that we could do variant calling for, named **vars_sra.txt**.

```bash

    while read i; do
        bash Vars_vcf_to_fa.sh $i
    done < Source/vars_sra.txt

```

## Formatting output data
With all variant calling finished, the individual DNA and peptide sequences were collected into FASTA and CSV files to be read later downstream.

```bash

    bash format_output.sh

```
