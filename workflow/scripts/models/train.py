"""
Generic script for training models
"""

method = snakemake.wildcards["method"]

train_data_file = snakemake.input["train"]

# TODO: train using method-specific function

model_file = snakemake.output["model"]
