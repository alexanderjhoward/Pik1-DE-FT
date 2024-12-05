# Specify GPU
import torch
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Set model and seed
model_checkpoint = "facebook/esm2_t6_8M_UR50D"
seed = 42

# Define averaged embedding function
def embeddings(dataloader):
    import numpy as np
    model.eval() # set model to evaluation mode to keep model parameters locked during evaluation
    embeddings = []

    for batch in dataloader: # loop over batches from val_dataloader
        batch = tuple(b.to(device) for b in batch) # move batch to GPU and create a dictionary of IDs, attention masks, and labels
        inputs = {"input_ids": batch[0],
                  "attention_mask": batch[1]
                  }
        with torch.no_grad():
            outputs = model(**inputs, output_hidden_states=True) # add inputs to model
        output_embeddings = outputs[-1][-1].detach().cpu().numpy().mean(axis=1) # move embeddings to CPU and append the mean along axis 1
        embeddings.append(output_embeddings)

    embeddings = np.concatenate(embeddings, axis=0) # reshape prediction and true_values to get single arrays
    return embeddings # return desired evaluation metrics

# Loop over Avr-PikC and Avr-PikF datasets to extract sequence embeddings
ligands = ["AvrPikC", "AvrPikF"]
for x in ligands:

    # Split data into training and validation sets
    import pandas as pd
    from sklearn.model_selection import train_test_split
    facs = pd.read_csv(f'../../1_Format_DE_Data/Output/Data/{x}_Enrichment.csv')
    train, val = train_test_split(facs, train_size=0.95, random_state=seed, stratify=facs['bin'])
    train = train[['hma', 'enrichment']]
    val = val[['hma', 'enrichment']]

    # Separate sequences and labels in each set
    train_seq = train['hma']
    val_seq = val['hma']

    # Import ESM tokenizer and model
    from transformers import AutoTokenizer, AutoModelForSequenceClassification
    tokenizer = AutoTokenizer.from_pretrained(model_checkpoint)
    model = AutoModelForSequenceClassification.from_pretrained(model_checkpoint, num_labels=1)
    model.to(device)

    # Tokenize and encode sequences
    max_seq_len = 80
    ## Training
    train_tokens = tokenizer.batch_encode_plus(
        train_seq.tolist(),
        max_length=max_seq_len,
        padding='longest',
        truncation=True,
        return_token_type_ids=False
    )
    train_seq = torch.tensor(train_tokens['input_ids'])
    train_mask = torch.tensor(train_tokens['attention_mask'])
    ## Validation
    val_tokens = tokenizer.batch_encode_plus(
        val_seq.tolist(),
        max_length=max_seq_len,
        padding='longest',
        truncation=True,
        return_token_type_ids=False
    )
    val_seq = torch.tensor(val_tokens['input_ids'])
    val_mask = torch.tensor(val_tokens['attention_mask'])

    # Define DataLoader
    from torch.utils.data import TensorDataset, DataLoader, SequentialSampler
    batch_size = 32
    ## Training
    train_data_wrap = TensorDataset(train_seq, train_mask)
    train_dataloader = DataLoader(train_data_wrap, sampler=SequentialSampler(train_data_wrap), batch_size=batch_size)
    ## Validation
    val_data_wrap = TensorDataset(val_seq, val_mask)
    val_dataloader = DataLoader(val_data_wrap, sampler=SequentialSampler(val_data_wrap), batch_size=batch_size)

    # Extract embeddings
    ## Training
    train_embeddings = embeddings(train_dataloader)
    train_output = pd.DataFrame.from_records(train_embeddings.reshape(train.shape[0], 320))
    train_output = pd.concat([train, train_output], axis=1)
    train_output.to_csv(f'../Output/{x}/train_embeddings.csv', index=False)
    ## Validation
    val_embeddings = embeddings(val_dataloader)
    val_output = pd.DataFrame.from_records(val_embeddings.reshape(val.shape[0], 320))
    val_output = pd.concat([val, val_output], axis=1)
    val_output.to_csv(f'../Output/{x}/val_embeddings.csv', index=False)