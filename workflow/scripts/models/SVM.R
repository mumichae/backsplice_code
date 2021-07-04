# Support Vector Machine on dataset of Wang et al. (2019), by the model of Wang et al. (2017)
library(e1071)
library(caret)
library(ggplot2)
library(pROC)
#library(DeepCirCode)

# Load datasets in lazy loading mode:
#data(HumanTrain)  # datasets are lazy loaded as variable "HumanTrain" until being used
#data(HumanTest)
#data(y_train_human)
#data(y_test_human)


freq_train <- readRDS(snakemake@input$train_features)
y_train <- data.frame(read.table(snakemake@input$train_labels, colClasses = 'character' ))
y_train <- cbind(as.integer(substr(y_train[, 1], 1, 1)), as.integer(substr(y_train[, 1], 2, 2)))

model_path <- snakemake@output$SVM_model

# we need to extract important sequence features! 
# frequencies of 1/2/3-mer compositions, normalized by length of their intron/exon -> 336 k-mer compositional features!!!

xy_train <- cbind(label = as.factor(y_train[,2]), data.frame(freq_train))


# Fit the RF model
set.seed(72)

svm_model <- svm(label ~., data = xy_train, method = "C-classification", kernel= "radial", gamma = 0.1, cost = 10) # gamma and cost may have to be tuned

summary(svm_model)

# svm_model$SV

#plot(svm_model, xy_train_human, most important vector, second most important vector, slice = list(?))

save(svm_model, file = model_path)

# Use model

#use fitted bagged model to predict label of new observation
#prediction_svm <- predict(svm_model, data.frame(freq_test))
#prediction_svm <- data.frame(cbind(label = y_test_human[, 2], prediction = prediction_svm, prediction_bin = round(prediction_svm, 0) ) ) 

#conf_table <- table(prediction_svm$label, prediction_svm$prediction_bin)
#conf_table


#roc <- roc(prediction_svm$label, prediction_svm$prediction)
#auc <- round(auc(prediction_svm$label, prediction_svm$prediction),4)

#create ROC plot
#plot <- ggroc(roc, colour = 'steelblue', size = 2) +
#  geom_abline(intercept = 1)+
#  ggtitle(paste0('SVM ROC Curve ', '(AUC = ', auc, ')'))+
#  theme_minimal()

#ggsave("ROC_SVM.png", plot)

# to load the model, use: 
# load("circRNA_RF.RData")



