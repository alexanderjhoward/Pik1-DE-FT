# Specify GPU
import torch
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Set model and seed to be used
model_checkpoint = "facebook/esm2_t6_8M_UR50D"

# Define model evaluation function
def evaluate(dataloader):
    import numpy as np
    model.eval() # set model to evaluation mode to keep model parameters locked during evaluation
    loss_val_total = 0 # initialize total validation loss, predicted logits, and true label variables
    predictions, true_vals = [], []

    for batch in dataloader: # loop over batches from val_dataloader
        batch = tuple(b.to(device) for b in batch) # move batch to GPU and create a dictionary of IDs, attention masks, and labels
        inputs = {"input_ids": batch[0],
                  "attention_mask": batch[1],
                  "labels": batch[2]
                  }
        with torch.no_grad():
            outputs = model(**inputs) # add inputs to model
        loss = outputs[0] # save loss
        logits = outputs[1] # save logits
        loss_val_total += loss.item() # accumulate loss in total validation loss variable
        logits = logits.detach().cpu().numpy() # move logits and labels to CPU and append to prediction and true label variables
        label_ids = inputs["labels"].cpu().numpy()
        predictions.append(logits)
        true_vals.append(label_ids)

    loss_val_avg = loss_val_total/len(dataloader) # once done with all batches, calcualte average validation loss
    predictions = np.concatenate(predictions, axis=0) # reshape prediction and true_values to get single arrays
    true_vals = np.concatenate(true_vals, axis=0)
    return loss_val_avg, predictions, true_vals # return desired evaluation metrics

# Loop over Avr-PikC and Avr-PikF datasets to generate predictions for both ligands
ligands = ["AvrPikC", "AvrPikF"]
for x in ligands:

    # Load in variant sequences
    import pandas as pd
    vars = pd.read_csv('../../3_Variant_Calling/Output/Pikp_HMA_Variants_Unique.csv')
    vars_seq = vars['SEQ'].tolist()
    vars['enrichment'] = 0 # placeholder enrichment score, true score is unknown
    vars_labels = vars['enrichment']

    # Define tokenizer to use
    from transformers import AutoTokenizer
    tokenizer = AutoTokenizer.from_pretrained(model_checkpoint)

    # Tokenize and format variant sequences for evaluation
    max_seq_len = 100
    tokens_vars = tokenizer.batch_encode_plus(
        vars_seq,
        max_length=max_seq_len,
        padding='longest',
        truncation=True,
        return_token_type_ids=False
    )
    vars_seq = torch.tensor(tokens_vars['input_ids'])
    vars_mask = torch.tensor(tokens_vars['attention_mask'])
    vars_y = torch.tensor(vars_labels.tolist())

    # Load into DataLoader
    from torch.utils.data import TensorDataset, DataLoader, RandomSampler, SequentialSampler
    batch_size = 32
    vars_data = TensorDataset(vars_seq, vars_mask, vars_y)
    vars_dataloader = DataLoader(vars_data, sampler=SequentialSampler(vars_data), batch_size=batch_size)

    # Apply finetuned weights to ESM-2
    from transformers import AutoModelForSequenceClassification
    model_save_path = f"../../2_Model_DE_Data/Output/{x}/esm2_t6_8M_UR50D_regression.model"
    model = AutoModelForSequenceClassification.from_pretrained(model_checkpoint, num_labels=1)
    model.load_state_dict(torch.load(model_save_path))
    model = model.to(device)

    # Evaluate variant data
    _, predictions, true_vals = evaluate(vars_dataloader)
    vars = pd.read_csv('../../3_Variant_Calling/Output/Pikp_HMA_Variants_Unique.csv')
    output_df = pd.concat([vars.reset_index(), pd.DataFrame(predictions, columns=[f'Predicted_{x}_Enrichment'])], axis=1, join='outer')
    output_df.to_csv(f'../Output/{x}_esm2_t6_8M_UR50D_predictions.csv', index=False)