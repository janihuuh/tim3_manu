## Init
import os
import numpy as np
import pandas as pd
import torch

from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.csv import CsvDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import VAE

# -------------------------------
# Parameters
# -------------------------------
use_cuda = True
n_epochs = 100
output_dir = "results/scvi/results"
input_dir = "results/scvi/input_files"

# Create results directory if missing
os.makedirs(output_dir, exist_ok=True)

# List of input CSVs (relative to input_dir)
input_files = [
    "files_here1.csv",
    "files_here2.csv"
]

# -------------------------------
# Load datasets
# -------------------------------
datasets = []
for fname in input_files:
    datasets.append(
        CsvDataset(
            filename=os.path.join(input_dir, fname),
            save_path="",
            sep=",",
            new_n_genes=False,
        )
    )

# Combine into one GeneExpressionDataset
all_dataset = GeneExpressionDataset()
all_dataset.populate_from_per_batch_list(Xs=[ds.X for ds in datasets])

# -------------------------------
# Train VAE
# -------------------------------
vae = VAE(
    all_dataset.nb_genes,
    n_batch=all_dataset.n_batches,
    n_labels=all_dataset.n_labels,
    n_hidden=128,
    n_latent=30,
    n_layers=2,
    dispersion="gene",
)

trainer = UnsupervisedTrainer(vae, all_dataset, train_size=1.0)
trainer.train(n_epochs=n_epochs)

# Save trained model
model_path = os.path.join(output_dir, "tim3_oneshot.pkl")
torch.save(trainer.model.state_dict(), model_path)

# -------------------------------
# Extract latent representation
# -------------------------------
full = trainer.create_posterior(trainer.model, all_dataset, indices=np.arange(len(all_dataset)))

latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

# Save embeddings
np.savetxt(os.path.join(output_dir, "tim3_oneshot_latent.csv"), latent, delimiter=",")
np.savetxt(os.path.join(output_dir, "tim3_oneshot_indices.csv"), batch_indices, delimiter=",")
