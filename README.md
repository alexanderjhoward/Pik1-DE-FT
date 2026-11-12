# Modeling *Pik-1* directed evolution data with fine-tuned ESM-2 models to phenotype naturally-evolved *Pik-1* variants
This repository contains all scripts used in the paper "Combining Directed Evolution with Machine Learning Enables Accurate Genotype-to-Phenotype Predictions".

## Data
Source data can be found at Zenodo (DOI:10.5281/zenodo.24174238). To reproduce our modeling/figures, download our data from Zenodo and place in the following directories:
* Pik1_DEFT_Source_Data.zip (unzipped) in Data/source/
* YSD_Analysis_Source_Data.zip (unzipped) in YSD_Analysis/Data/

## Analysis
All analyses were completed in R or Python. RStudio was used to run the R scripts and Google Collab (T4 GPU) was used to run the Python scripts.

## Outputs
Model embeddings, weights, and figures can be found at Zenodo (DOI:10.5281/zenodo.24174238). To collect these script outputs, download our data from Zenodo and place in the following directories:
* Pik1_DEFT_Outputs.zip (unzipped) in Outputs/

## Overview
1. [Pik-1 variant calling](https://github.com/alexanderjhoward/Pik1-DE-FT/tree/main/3kRGP): The 3,000 Rice Genomes Project dataset was searched for uncharacterized variants of the Pik-1 HMA domain to phenotype with our modeling.
2. [Model directed evolution data](https://github.com/alexanderjhoward/Pik1-DE-FT/tree/main/Pik1_DE-FT.ipynb): The protein language model ESM-2 was used to model Pik-1 binding against Avr-PikC and Avr-PikF (as measured through prior directed evolution of the Pik-1 receptor). This modeling approach was also tested on the human enzyme NUDT15. Model approaches are compared and visualized with this script.
3. [YSD Pik-1 phenotyping](https://github.com/alexanderjhoward/Pik1-DE-FT/tree/main/YSD_Analysis): Pik-1 HMA domains of interest were tested for binding against Avr-Pik in YSD and imaged/sorted to quantify relative binding against each effector variant. Phenotyping results are visualized with this script.
