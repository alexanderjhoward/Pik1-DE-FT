# Recapitulating genotype-to-phenotype approach by modeling NUDT15 mutagenesis data

## Background
In order to test the generalizability of our modeling approach for genotype-to-phenotype predictions, we wanted to model an additional mutagenesis dataset which connects protein variation to phenotypic effects. For this project, we chose the [Suiter *et al.*, 2020](https://www.pnas.org/doi/10.1073/pnas.1915680117) paper exploring the effect of NUDT15 variation on enzyme stability and functionality. This dataset was interesting because loss-of-function NUDT15 variants have been previously reported to cause cytotoxicity in patients treated with thiopurine drugs. We wanted to explore how well we could correlate NUDT15 sequence varation with enzyme functionality (and in extension, cytotoxicity risk) by fine-tuning ESM-2. This model would then be used to characterize known clinical variants of the enzyme and any currently uncharacterized variants.

## Method
The **nudt15.R** script was used to generate a list of clinical variants and uncharacterized genomAD variants of NUDT15 for downstream phenotyping, as well as process the rest of the Suiter *et al.* data for modeling.

The **esm2_nudt15_finetuning.py** script was used to generate fine-tuned ESM-2 models from the NUDT15 assay functional scores.
