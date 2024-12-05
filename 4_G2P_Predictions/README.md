# Utilizing fine-tuned ESM-2 models to conduct genotype-to-phenotype analyses on Avr-Pik ligand binding 

## Background
After we identifed uncharacterized *Pik-1* ligand binding domains, we wanted to see what our fine-tuned models would predict for these variants in terms of Avr-PikC/Avr-PikF ligand binding. These predictions could in turn inform the selection of variants of interest to test more in-depth downstream.  

## Method
The **esm2_g2p_predictions.py** script was used to predict the Avr-PikC/Avr-PikF pES values of our 3k RGP *Pik-1* variants. 
