# Random Forest on dataset of Wang et al. (2019), by the model of Wang et al. (2017)
library(randomForest)
library(caret)
#library(DeepCirCode)
library(ggplot2)
library(pROC)

freq_test <- readRDS(snakemake@input$test_features)
y_test <- data.frame(read.table(snakemake@input$test_labels, colClasses = 'character' ))
#y_train <- data.frame(read.table('//wsl$/Ubuntu-20.04/home/laura/ML4RG_project/backsplice_code/results/processed_data/features/DeepCirCode/Wang2019/y_matrix.txt',
#                                 colClasses = 'character'))
y_test <- cbind(as.integer(substr(y_test[, 1], 1, 1)), as.integer(substr(y_test[, 1], 2, 2)))

rf_model <- get(load(snakemake@input$RF_model))
print(rf_model)
plot_path <- snakemake@output$plot
prediction_path <- snakemake@output$prediction


# Use model

#use fitted bagged model to predict value of new observation
prediction_rf <- predict(rf_model, newdata = data.frame(freq_test))
prediction_rf <- data.frame(cbind(label = y_test[, 2], prediction = prediction_rf, prediction_bin = round(prediction_rf, digits = 0)))

conf_table <- table(prediction_rf$label, prediction_rf$prediction_bin)

conf_table

write.table(prediction_rf, file = prediction_path)

roc <- roc(prediction_rf$label, prediction_rf$prediction)
auc <- round(auc(prediction_rf$label, prediction_rf$prediction), 4)

#create ROC plot
plot <- ggroc(roc, colour = 'steelblue', size = 2) +
  geom_abline(intercept = 1) +
  ggtitle(paste0('RF ROC Curve ', '(AUC = ', auc, ')')) +
  theme(plot.background = element_rect(fill = "white"))


ggsave(plot_path, plot)


