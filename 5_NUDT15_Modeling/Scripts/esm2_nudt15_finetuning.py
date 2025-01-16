# Specify GPU
import torch
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Set model and seed to be used
model_checkpoint = "facebook/esm2_t6_8M_UR50D"
seed = 42

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

# Split data into training and validation sets
import pandas as pd
from sklearn.model_selection import train_test_split
facs = pd.read_csv('../Output/Data/NUDT15_cleaned.csv')
train, val = train_test_split(facs, train_size=0.90, random_state=seed, stratify=facs['bin'])

# Separate sequences and labels
train_seq = train['Protein'].tolist()
val_seq = val['Protein'].tolist()
train_labels = train['Score']
val_labels = val['Score']

# Tokenize sequences
from transformers import AutoTokenizer
tokenizer = AutoTokenizer.from_pretrained(model_checkpoint)
train_tokenized = tokenizer(train_seq)
val_tokenized = tokenizer(val_seq)

# Load into datasets
from datasets import Dataset
train_dataset = Dataset.from_dict(train_tokenized)
val_dataset = Dataset.from_dict(val_tokenized)
train_dataset = train_dataset.add_column("labels", train_labels)
val_dataset = val_dataset.add_column("labels", val_labels)

# Define evaluation metric
from evaluate import load
metric = load('spearmanr')
def compute_metrics(eval_pred):
    predictions, labels = eval_pred
    return metric.compute(predictions=predictions, references=labels)

# Load model
from transformers import AutoModelForSequenceClassification, TrainingArguments, Trainer, set_seed
import math
set_seed(seed)
model = AutoModelForSequenceClassification.from_pretrained(model_checkpoint, num_labels=1)  # Doing regression, only need one output
model_name = model_checkpoint.split("/")[-1]
batch_size = math.floor(train.shape[0] / 500)
args = TrainingArguments(
    f"{model_name}-finetuned-regression",
    evaluation_strategy="epoch",
    save_strategy="epoch",
    learning_rate=1e-05,
    per_device_train_batch_size=batch_size,
    per_device_eval_batch_size=batch_size,
    num_train_epochs=20,
    weight_decay=0.01,
    seed=seed,
    load_best_model_at_end=True,
    metric_for_best_model="spearmanr",
)
trainer = Trainer(
    model,
    args,
    train_dataset=train_dataset,
    eval_dataset=val_dataset,
    tokenizer=tokenizer,
    compute_metrics=compute_metrics
)

# Train model
trainer.train()
model_save_path = "../Output/esm2_t6_8M_UR50D_regression.model"
torch.save(model.state_dict(), model_save_path)

# Import variant data to phenotype
vars = pd.read_csv('../Output/Data/NUDT15_variants.csv')
vars['Score'] = 0
vars_seq = vars['seq'].tolist()
vars_labels = vars['Score']

# Tokenize and format validation sequences for evaluation
max_seq_len = 180  
tokens_val = tokenizer.batch_encode_plus(
    val_seq,
    max_length=max_seq_len,
    padding='longest',
    truncation=True,
    return_token_type_ids=False
)
tokens_vars = tokenizer.batch_encode_plus(
    vars_seq,
    max_length = max_seq_len,
    padding='longest',
    truncation=True,
    return_token_type_ids=False
)
val_seq = torch.tensor(tokens_val['input_ids'])
val_mask = torch.tensor(tokens_val['attention_mask'])
val_y = torch.tensor(val_labels.tolist())
vars_seq = torch.tensor(tokens_vars['input_ids'])
vars_mask = torch.tensor(tokens_vars['attention_mask'])
vars_y = torch.tensor(vars_labels.tolist())

# Load into DataLoader
from torch.utils.data import TensorDataset, DataLoader, RandomSampler, SequentialSampler
batch_size = 32
val_data = TensorDataset(val_seq, val_mask, val_y)
val_dataloader = DataLoader(val_data, sampler=SequentialSampler(val_data), batch_size=batch_size)
vars_data = TensorDataset(vars_seq, vars_mask, vars_y)
vars_dataloader = DataLoader(vars_data, sampler=SequentialSampler(vars_data), batch_size=batch_size)

# Apply finetuned weights to ESM-2
model = AutoModelForSequenceClassification.from_pretrained(model_checkpoint, num_labels=1)
model.load_state_dict(torch.load(model_save_path))
model = model.to(device)

# Evaluate samples on validation data
_, predictions, true_vals = evaluate(val_dataloader)
_, vars_predictions, vars_true_vals = evaluate(vars_dataloader)
output_df = pd.concat([val.reset_index(), pd.DataFrame(predictions, columns=['Predicted_Score'])], axis=1, join='outer')
output_df.to_csv('../Output/esm2_t6_8M_UR50D_regression_predictions.csv', index=False)
vars_output_df = pd.concat([vars.reset_index(), pd.DataFrame(vars_predictions, columns = ['Predicted_Score'])], axis=1, join='outer')
vars_output_df.to_csv('../Output/esm2_t6_8M_UR50D_regression_variant_predictions.csv', index=False)
