# Title     : Train model
# Objective : Generic script for training models in R
# Created by: mumichae
# Created on: 7/6/21

method <- snakemake@params$method
train_data_file <- snakemake@input$train

model_file <- snakemake@output$model
