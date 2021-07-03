# Random Forest on dataset of Wang et al. (2019), by the model of Wang et al. (2017)
library(randomForest)
library(caret)
#library(DeepCirCode)
library(ggplot2)
library(pROC)

freq_train <- readRDS(snakemake@input$train_features)
y_train <- data.frame(read.table(snakemake@input$train_labels, colClasses = 'character' ))
#y_train <- data.frame(read.table('//wsl$/Ubuntu-20.04/home/laura/ML4RG_project/backsplice_code/results/processed_data/features/DeepCirCode/Wang2019/y_matrix.txt',
#                                 colClasses = 'character'))
y_train <- cbind(as.integer(substr(y_train[, 1], 1, 1)), as.integer(substr(y_train[, 1], 2, 2)))

model_path <- snakemake@output$RF_model
#freq_test <- readRDS(snakemake@input$test)
#plot_path <- snakemake@output$plot

# we need to extract important sequence features! 
# frequencies of 1/2/3-mer compositions, normalized by length of their intron/exon -> 336 k-mer compositional features!!!
xy_train <- cbind(label = y_train[, 2], data.frame(freq_train))

# Fit the RF model
set.seed(72)
rf_model <- randomForest(
  formula = label ~ .,
  data = xy_train
)
# display model
rf_model

#find Number of Trees that produce lowest test MSE
which.min(rf_model$mse)

# find RMSE of best model, here: average difference between predicted label and observed label
sqrt(rf_model$mse[which.min(rf_model$mse)])


# plot test MSE by number of trees
plot(rf_model)

# plot importance of each predictor variable
varImpPlot(rf_model)


save(rf_model, file = model_path)

# Use model

# #use fitted bagged model to predict Ozone value of new observation
# prediction_rf <- predict(rf_model, newdata = data.frame(freq_test))
# prediction_rf <- data.frame(cbind(label = y_test_human[, 2], prediction = prediction_rf, prediction_bin = round(prediction_rf, digits = 0)))
# 
# conf_table <- table(prediction_rf$label, prediction_rf$prediction_bin)
# 
# conf_table
# 
# 
# roc <- roc(prediction_rf$label, prediction_rf$prediction)
# auc <- round(auc(prediction_rf$label, prediction_rf$prediction), 4)
# 
# #create ROC plot
# plot <- ggroc(roc, colour = 'steelblue', size = 2) +
#   geom_abline(intercept = 1) +
#   ggtitle(paste0('RF ROC Curve ', '(AUC = ', auc, ')')) +
#   theme_minimal()
# 
# ggsave(plot_path, plot)


