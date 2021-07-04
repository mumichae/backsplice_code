# assess performance of our models
library(ggplot2)
library(pROC)
library(data.table)

#read in prediction
RF_p <- data.table(read.table(snakemake@input$RF_prediction, header=TRUE))
SVM_p <-  data.table(read.table(snakemake@input$SVM_prediction, header = TRUE))
#DCC_p <- read.table(snakemake@input$DCC_prediction)
#JEDI_p <-  read.table(snakemake@input$JEDI_prediction)


# read output paths
barplot_path <- snakemake@output$barplot
roc_path <- snakemake@output$roc
#prc_path <- snakemake@output$prc

# confusion matrices
RF_conf <- table(RF_p$label, RF_p$prediction_bin)
RF_conf
SVM_conf <- table(SVM_p$label, SVM_p$prediction_bin)
SVM_conf


general_performance <- data.table(model = c("RandomForest", "SVM")) #, "Random", "DCC", "JEDI"))
# acc = TP / (TP+FP)
acc <- c(RF_conf[2,2] / sum(RF_conf[,2]), 
         SVM_conf[2,2] / sum(SVM_conf[,2]))
# balanced acc = 
# sens = TP / (TP+FN)
sens <- c(RF_conf[2,2] / sum(RF_conf[2,]), 
         SVM_conf[2,2] / sum(SVM_conf[2,]))
# F1 = 2*Sens*Acc / (Sens+Acc)
f1 <- (2*sens*acc)/(sens+acc)
# MCC = (TP*TN - FP*FN) / (sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN)))
# mcc = c((RF_conf[2,2]*RF_conf[1,1] - RF_conf[1,2]*RF_conf[2,1]) / (sqrt(sum(RF_conf[, 2]) * sum(RF_conf[2,]) * sum(RF_conf[1,]) * sum(RF_conf[,1]))) ,
#        (SVM_conf[2,2]*SVM_conf[1,1] - SVM_conf[1,2]*SVM_conf[2,1]) / (sqrt(sum(SVM_conf[, 2]) * sum(SVM_conf[2,]) * sum(SVM_conf[1,]) * sum(SVM_conf[,1]))))


general_performance <- cbind(general_performance, Accuracy = acc, Sensitivity = sens, 
                             F1 = f1)
general_performance

general_performance <- melt(general_performance, id.vars = "model")
general_performance


# accuracy/Q2, sensitivity, MCC, F1
# precision recall curves, roc curves
plot <- ggplot(general_performance)+
  geom_bar(mapping = aes(x = variable, y = value, fill = model), stat = "identity", position = position_dodge2())+
  geom_text(aes(x = variable, y = -0.02, label = round(value, 2)),
            position = position_dodge2(width = 0.9), size = 3.5)+
  ggtitle("Performance Assessment")+
  xlab("metric")+
  ylab("performance")+
  #scale_fill_manual(name="Model", values=c("#66FF66", "#FF6600", "#CC0000", "#CC6699", "#9999CC"))+
  theme_bw()

barplot_path
ggsave( barplot_path, plot)


#roc curve
roc <- roc(RF_p$label, RF_p$prediction)
auc <- round(auc(RF_p$label, RF_p$prediction), 4)

#create ROC plot
roc_plot <- ggroc(roc, colour = 'steelblue', size = 2) +
  geom_abline(intercept = 1) +
  ggtitle(paste0('RF ROC Curve ', '(AUC = ', auc, ')')) +
  theme(plot.background = element_rect(fill = "white"))


ggsave(roc_path, roc_plot)

#TODO: precision-recall curves

