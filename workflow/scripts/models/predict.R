# Title     : Predict on test data
# Objective : Generic script for predicting models in R
# Created by: mumichae
# Created on: 7/6/21

method <- snakemake@wildcards$method

model_file <- snakemake@input$model
test_data_file <- snakemake@input$test

# TODO: predict R models

prediction_file <- snakemake@output$prediction
