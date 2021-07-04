# Support Vector Machine on dataset of Wang et al. (2019), by the model of Wang et al. (2017)
library(e1071)
library(caret)
library(ggplot2)
library(pROC)
#library(DeepCirCode)

freq_test <- readRDS(snakemake@input$test_features)
y_test <- data.frame(read.table(snakemake@input$test_labels, colClasses = 'character' ))
#y_train <- data.frame(read.table('//wsl$/Ubuntu-20.04/home/laura/ML4RG_project/backsplice_code/results/processed_data/features/DeepCirCode/Wang2019/y_matrix.txt',
#                                 colClasses = 'character'))
y_test <- cbind(as.integer(substr(y_test[, 1], 1, 1)), as.integer(substr(y_test[, 1], 2, 2)))

svm_model <- get(load(snakemake@input$model))

plot_path <- snakemake@output$plot
prediction_path <- snakemake@output$prediction


# Use model

#use fitted bagged model to predict label of new observation
prediction_svm <- predict(svm_model, data.frame(freq_test))
prediction_svm <- data.frame(cbind(label = y_test[, 2], prediction = prediction_svm, prediction_bin = round(as.numeric(prediction_svm), 0) ) ) 

conf_table <- table(prediction_svm$label, prediction_svm$prediction)
conf_table

write.table(prediction_svm, file = prediction_path)

roc <- roc(prediction_svm$label, prediction_svm$prediction)
auc <- round(auc(prediction_svm$label, prediction_svm$prediction),4)

#create ROC plot
plot <- ggroc(roc, colour = 'steelblue', size = 2) +
  geom_abline(intercept = 1)+
  ggtitle(paste0('SVM ROC Curve ', '(AUC = ', auc, ')')) +
  theme(plot.background = element_rect(fill = "white"))

ggsave(plot_path, plot)




