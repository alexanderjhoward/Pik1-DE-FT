# Modeling *Pik-1* directed evolution data with fine-tuned ESM-2 models to phenotype naturally-evolved *Pik-1* variants
This repository contains all scripts and data necessary to reproduce the analysis conducted in the paper "Combining Directed Evolution with Machine Learning Enables Accurate Genotype-to-Phenotype Predictions on Rare Receptor-Ligand Recognition Profiles".

## Overview
1. [Format directed evolution data](https://github.com/alexanderjhoward/Pik1-DE-FT/tree/main/1_Format_DE_Data): Sequencing data from our starting yeast surface display library and post-selection libraries was cleaned and processed to calculate an enrichment score (ES) for each *Pik-1* sequence variant which quantifies the relative change in variant abundance following fluorescence-activated cell sorting.
2. [Model directed evolution data](https://github.com/alexanderjhoward/Pik1-DE-FT/tree/main/2_Model_DE_Data): The protein language model ESM-2 was fine-tuned on our *Pik-1* directed evolution data to calculate a predicted enrichment score (pES) for input *Pik-1* sequence variants.
3. [Variant calling](https://github.com/alexanderjhoward/Pik1-DE-FT/tree/main/3_Variant_Calling): The 3,000 Rice Genomes Project dataset was searched for uncharacterized variants of the *Pik-1* ligand binding domain to phenotype with our modeling.
4. [Genotype-to-phenotype predictions](https://github.com/alexanderjhoward/Pik1-DE-FT/tree/main/4_G2P_Predictions): Our fine-tuned ESM-2 models were used to calculate pES values for the naturally-evolved variants of *Pik-1*.
5. [Model NUDT15 data](https://github.com/alexanderjhoward/Pik1-DE-FT/tree/main/5_NUDT15_Modeling): Recapitulate our genotype-to-phenotype approach with another dataset (NUDT15 enzyme functionality).
6. [Visualize results](https://github.com/alexanderjhoward/Pik1-DE-FT/tree/main/6_Visualize_Results): The modeling and variant calling results were visualized.
