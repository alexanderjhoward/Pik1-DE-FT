# Cleaning Pikh-1 directed evolution data for downstream modeling

## Background
After applying directed evolution to the *Pikh-1* ligand binding domain, next-generation sequencing was used on our YSD libraries to survey the change in sequence abundance before and after FACS selection. Here we process that data to obtain an "enrichment score" that quantifies variant fitness towards Avr-PikC and Avr-PikF ligand binding. 

## Method
I ran the **clean_facs_data.R** script to clean and score sequence variants for downstream modeling.
