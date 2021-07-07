suppressPackageStartupMessages({
  library(keras)
  library(reticulate)
  library(DeepCirCode)
  library(data.table)
})

# check of keras package
# model <- keras_model_sequential()

model <- load_model_hdf5(filepath = snakemake@input$model)
# model_path <- system.file("DeepCirCode_bestmodel_human.hdf5", package ="DeepCirCode")
# model <- load_model_hdf5(filepath = model_path)

test_set <- fread(snakemake@input$test_features)

##Converting one-hot encoded sequence into 3D matrix as the input for DeepCirCode
x_test <- convStringToMatrix(test_set$encoded_seq)
y_test <- to_categorical(test_set$label, 2)

##Testing DeepCirCode with test sets
# preds <- DeepCirCode_test(x_test, y_test, "human")
pred_df <- data.table(
  label = y_test[, 2],
  score = model %>% predict_proba(x_test),
  prediction = model %>% predict_classes(x_test)
)

fwrite(pred_df, snakemake@output$prediction, sep = '\t')
