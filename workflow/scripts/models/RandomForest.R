# Random Forest on dataset of Wang et al. (2019), by the model of Wang et al. (2017)
library(randomForest)
library(caret)
library(DeepCirCode)
library(ggplot2)
library(pROC)

freq_train <- readRDS(snakemake@input$train)

# we need to extract important sequence features! 
# frequencies of 1/2/3-mer compositions, normalized by length of their intron/exon -> 336 k-mer compositional features!!!
xy_train_human <- cbind(label = y_train_human[, 2], data.frame(freq_train))

# Fit the RF model
set.seed(72)
rf_model <- randomForest(
  formula = label ~ .,
  data = xy_train_human
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


save(rf_model, file = "circRNA_RF.RData")

# Use model

#use fitted bagged model to predict Ozone value of new observation
prediction_rf <- predict(rf_model, newdata = data.frame(freq_test))
prediction_rf <- data.frame(cbind(label = y_test_human[, 2], prediction = prediction_rf, prediction_bin = round(prediction_rf, digits = 0)))

conf_table <- table(prediction_rf$label, prediction_rf$prediction_bin)

conf_table


roc <- roc(prediction_rf$label, prediction_rf$prediction)
auc <- round(auc(prediction_rf$label, prediction_rf$prediction), 4)

#create ROC plot
plot <- ggroc(roc, colour = 'steelblue', size = 2) +
  geom_abline(intercept = 1) +
  ggtitle(paste0('RF ROC Curve ', '(AUC = ', auc, ')')) +
  theme_minimal()

ggsave("ROC_RF.png", plot)


