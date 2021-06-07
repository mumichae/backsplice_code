"""
Generic script for predicting models
"""

method = snakemake.wildcards['method']

model_file = snakemake.input['model']
test_data_file = snakemake.input['test']

# TODO: predict specifically for method

prediction_file = snakemake.output['prediction']
