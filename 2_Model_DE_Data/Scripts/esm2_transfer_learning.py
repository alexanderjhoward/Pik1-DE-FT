# Specify GPU
import torch
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Set model and seed to be used
model_checkpoint = "facebook/esm2_t6_8M_UR50D"
seed = 42

# Create an empty dataframe to save results to
import pandas as pd
df = pd.DataFrame(columns=['Model', 'Ligand', 'Spearman', 'Significance'])
df2 = pd.DataFrame(columns=['Model', 'Ligand', 'Spearman', 'Significance'])

# Loop over Avr-PikC and Avr-PikF datasets to generate fine-tuned ESM-2 models
ligands = ["AvrPikC", "AvrPikF"]
for x in ligands:

    # Import embedding data
    train = pd.read_csv(f'../Output/{x}/train_embeddings.csv')
    val = pd.read_csv(f'../Output/{x}/val_embeddings.csv')
    train_cs = pd.read_csv(f'../Output/{x}/train_embeddings_cluster_split.csv')
    val_cs = pd.read_csv(f'../Output/{x}/val_embeddings_cluster_split.csv')
    
    # Separate features and labels
    X_train = train.drop(['hma', 'enrichment'], axis=1)
    y_train = train['enrichment']
    X_validation = val.drop(['hma', 'enrichment'], axis=1)
    y_validation = val['enrichment']
    X_train_cs = train_cs.drop(['hma', 'enrichment'], axis=1)
    y_train_cs = train_cs['enrichment']
    X_validation_cs = val_cs.drop(['hma', 'enrichment'], axis=1)
    y_validation_cs = val_cs['enrichment']

    ### Gradient-boosted decision tree
    # Define model
    from catboost import CatBoostRegressor, Pool, metrics, cv
    from sklearn.metrics import accuracy_score
    params = {
        'iterations': 1000,
        'random_seed': seed,
        'use_best_model': True,
        'od_type': 'Iter',
        'od_wait': 40,  # early stop after 40 rounds of non-improvement of validation loss
        'task_type': 'GPU'
    }
    train_pool = Pool(X_train, y_train)
    validate_pool = Pool(X_validation, y_validation)
    train_cs_pool = Pool(X_train_cs, y_train_cs)
    validate_cs_pool = Pool(X_validation_cs, y_validation_cs)
    model = CatBoostRegressor(**params)

    ## Random split
    # Fit model
    model.fit(train_pool, eval_set=validate_pool)

    # Predict with model
    predictions = model.predict(X_validation)
    cb_output_df = pd.concat([val.reset_index(), pd.DataFrame(predictions, columns=['Predicted_Enrichment'])], axis=1, join='outer')

    # Calculate Spearman R
    from scipy.stats import spearmanr
    cb_r,cb_p = spearmanr(cb_output_df['enrichment'], cb_output_df['Predicted_Enrichment'])
    
    # Calculate RMSD
    from sklearn.metrics import root_mean_squared_error
    cb_rmse = root_mean_squared_error(cb_output_df['enrichment'], cb_output_df['Predicted_Enrichment'])
    df = pd.concat([df, pd.DataFrame({'Model': "CatBoost", 'Ligand': x, 'Spearman': [cb_r], 'Significance': [cb_p], 'RMSE': [cb_rmse]})], ignore_index=True)

    ## Cluster split
    # Fit model
    model.fit(train_cs_pool, eval_set=validate_cs_pool)

    # Predict with model
    predictions = model.predict(X_validation_cs)
    cb_output_df = pd.concat([val.reset_index(), pd.DataFrame(predictions, columns=['Predicted_Enrichment'])], axis=1, join='outer')

    # Calculate Spearman R
    from scipy.stats import spearmanr
    cb_r,cb_p = spearmanr(cb_output_df['enrichment'], cb_output_df['Predicted_Enrichment'])
    
    # Calculate RMSD
    from sklearn.metrics import root_mean_squared_error
    cb_rmse = root_mean_squared_error(cb_output_df['enrichment'], cb_output_df['Predicted_Enrichment'])
    df2 = pd.concat([df2, pd.DataFrame({'Model': "CatBoost", 'Ligand': x, 'Spearman': [cb_r], 'Significance': [cb_p], 'RMSE': [cb_rmse]})], ignore_index=True)
    
    ### Linear regression
    # Define model
    from sklearn.linear_model import ElasticNetCV
    model = ElasticNetCV(random_state=seed)

    ## Random split
    # Fit model
    model.fit(X_train, y_train)

    # Predict with model
    predictions = model.predict(X_validation)
    lr_output_df = pd.concat([val.reset_index(), pd.DataFrame(predictions, columns=['Predicted_Enrichment'])], axis=1, join='outer')

    # Calculate Spearman R
    lr_r,lr_p = spearmanr(lr_output_df['enrichment'], lr_output_df['Predicted_Enrichment'])
    
    # Calculate RMSD
    lr_rmse = root_mean_squared_error(lr_output_df['enrichment'], lr_output_df['Predicted_Enrichment'])
    df = pd.concat([df, pd.DataFrame({'Model': "ElasticNet Regression", 'Ligand': x, 'Spearman': [lr_r], 'Significance': [lr_p], 'RMSE': [lr_rmse]})], ignore_index=True)

    ## Cluster split
    # Fit model
    model.fit(X_train_cs, y_train_cs)

    # Predict with model
    predictions = model.predict(X_validation_cs)
    lr_output_df = pd.concat([val.reset_index(), pd.DataFrame(predictions, columns=['Predicted_Enrichment'])], axis=1, join='outer')

    # Calculate Spearman R
    lr_r,lr_p = spearmanr(lr_output_df['enrichment'], lr_output_df['Predicted_Enrichment'])
    
    # Calculate RMSD
    lr_rmse = root_mean_squared_error(lr_output_df['enrichment'], lr_output_df['Predicted_Enrichment'])
    df2 = pd.concat([df2, pd.DataFrame({'Model': "ElasticNet Regression", 'Ligand': x, 'Spearman': [lr_r], 'Significance': [lr_p], 'RMSE': [lr_rmse]})], ignore_index=True)
    
    ### SVR
    # Define model
    from sklearn import svm
    model = svm.SVR()

    ## Random split
    # Fit model
    model.fit(X_train, y_train)

    # Predict with model
    predictions = model.predict(X_validation)
    svr_output_df = pd.concat([val.reset_index(), pd.DataFrame(predictions, columns=['Predicted_Enrichment'])], axis=1, join='outer')

    # Calculate Spearman R
    svr_r,svr_p = spearmanr(svr_output_df['enrichment'], svr_output_df['Predicted_Enrichment'])
    
    # Calculate RMSD
    svr_rmse = root_mean_squared_error(svr_output_df['enrichment'], svr_output_df['Predicted_Enrichment'])
    df = pd.concat([df, pd.DataFrame({'Model': "SVR", 'Ligand': x, 'Spearman': [svr_r], 'Significance': [svr_p], 'RMSE': [svr_rmse]})], ignore_index=True)

    ## Cluster split
    # Fit model
    model.fit(X_train_cs, y_train_cs)

    # Predict with model
    predictions = model.predict(X_validation_cs)
    svr_output_df = pd.concat([val.reset_index(), pd.DataFrame(predictions, columns=['Predicted_Enrichment'])], axis=1, join='outer')

    # Calculate Spearman R
    svr_r,svr_p = spearmanr(svr_output_df['enrichment'], svr_output_df['Predicted_Enrichment'])
    
    # Calculate RMSD
    svr_rmse = root_mean_squared_error(svr_output_df['enrichment'], svr_output_df['Predicted_Enrichment'])
    df2 = pd.concat([df2, pd.DataFrame({'Model': "SVR", 'Ligand': x, 'Spearman': [svr_r], 'Significance': [svr_p], 'RMSE': [svr_rmse]})], ignore_index=True)
    
    ### ESM-2 fine-tuned model
    ## Random split
    # Load in predictions
    ft_output_df = pd.read_csv(f'../Output/{x}/esm2_t6_8M_UR50D_regression_predictions.csv')

    # Calculate Spearman R
    ft_r,ft_p = spearmanr(ft_output_df['enrichment'], ft_output_df['Predicted_Enrichment'])

    # Calculate RMSD
    ft_rmse = root_mean_squared_error(ft_output_df['enrichment'], ft_output_df['Predicted_Enrichment'])
    df = pd.concat([df, pd.DataFrame({'Model': "Finetuned ESM-2", 'Ligand': x, 'Spearman': [ft_r], 'Significance': [ft_p], 'RMSE': [ft_rmse]})], ignore_index=True)

    ## Cluster split
    # Load in predictions
    ft_output_df = pd.read_csv(f'../Output/{x}/esm2_t6_8M_UR50D_regression_predictions_cs.csv')

    # Calculate Spearman R
    ft_r,ft_p = spearmanr(ft_output_df['enrichment'], ft_output_df['Predicted_Enrichment'])

    # Calculate RMSD
    ft_rmse = root_mean_squared_error(ft_output_df['enrichment'], ft_output_df['Predicted_Enrichment'])
    df2 = pd.concat([df2, pd.DataFrame({'Model': "Finetuned ESM-2", 'Ligand': x, 'Spearman': [ft_r], 'Significance': [ft_p], 'RMSE': [ft_rmse]})], ignore_index=True)

# Save results
df.to_csv('../Output/table1_model_comparisons.csv', index=False)
df2.to_csv('../Output/table2_model_comparisons.csv', index=False)
