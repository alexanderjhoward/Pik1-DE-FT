# Modeling Pikh-1 directed evolution data

## Background
We wanted to explore how well we could correlate sequence varation with ligand binding through protein language modeling. For this project, we utilized the ESM-2 protein language model and compared the performance of models which were generated through fine-tuning versus transfer learning.  

## Method
The **esm2_average_feature_extraction.py** script was used to extract the embeddings of each sequence in our Avr-PikC and Avr-PikF data (maintaining the same train/validation split as our fine-tuned models for ease of cross comparison). These embeddings were averaged across the all amino acid positions to reduce the feature dimensionality.

The **cluster_split.R** script was used to manually do a train/validation split of the Avr-PikC and Avr-PikF data based on the embedding latent space.

The **esm2_finetuning.py** script was used to generate fine-tuned ESM-2 models from our Avr-PikC and Avr-PikF selection data on both a random train/validation split and the latent space train/validation split.

The **esm2_transfer_learning.py** script was used to conduct transfer learning with the previously generated embeddings, making a CatBoost model, an ElasticNet linear regression model, and a SVR model. This was done for both the same random train train/validation split used for the ESM-2 fine-tuning and the latent space train/validation split. The Spearman R correlation coefficient and RMSE for the validation predictions made by these models (as well as our fine-tuned models) was also calculated to compare each model's performance on the prediction task. 
