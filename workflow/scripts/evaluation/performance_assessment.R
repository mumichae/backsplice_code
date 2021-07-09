saveRDS(snakemake, ".snakemake/snakemake_eval.RDS")

# assess performance of our models
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(yardstick)
  library(dplyr)
  library(mltools)
})

# global variables
CONF_MAT_NAMES <- c('TP', 'FP', 'FN', 'TN')

pred_files <- snakemake@input$predictions
methods <- snakemake@params$method
sources <- snakemake@params$source

metrics_path <- snakemake@output$metrics
barplot_path <- snakemake@output$barplot
roc_path <- snakemake@output$roc
pr_path <- snakemake@output$pr

#read in prediction
predictions <- mapply(
  function(file, method_name, source_name) {
    dt <- fread(file)
    dt[, method := method_name]
    dt[, source := source_name]
    dt[, label := as.factor(label)]
    dt[, prediction := as.factor(prediction)]
    dt
  }, pred_files, methods, sources, SIMPLIFY = FALSE
)
print(paste('read', length(predictions), 'files'))

# confusion matrices
conf <- sapply(predictions, function(x) {
  table(x$label, x$prediction)
})
conf <- t(conf)
colnames(conf) <- CONF_MAT_NAMES
rownames(conf) <- paste(methods, sources, sep = '_')

print('confidence matrix')
print(conf)

general_performance <- as.data.table(conf, keep.rownames = 'id')
general_performance[, source := sources]
general_performance[, method := methods]
# acc = (TP+TN)/(TP+TN+FP+FN)
general_performance[, Accuracy := (TP + TN) / (TP + TN + FP + FN)]
# balanced acc = (TP/(TP+FP) + TN/(TN+FN)) / 2
general_performance[, Balanced_Accuracy := ((TP/(TP+FP)) + (TN/(TN+FN))) / 2]
# spec = TP / (TP+FP)
general_performance[, Specificity := TP / (TP + FP)]
# sens = TP / (TP+FN)
general_performance[, Sensitivity := TP / (TP + FN)]
# F1 = 2*Sens*Spec / (Sens+Spec)
general_performance[, F1 := (2 * Sensitivity * Specificity) / (Sensitivity + Specificity)]
# MCC = (TP*TN - FP*FN) / (sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN)))
general_performance[, MCC := mapply(mcc, TP = TP, FP = FP, TN = TN, FN = FN)]
setorder(general_performance, id)

print('general performance')
print(head(general_performance))

fwrite(general_performance, metrics_path, sep = '\t')

per_melt <- melt(general_performance, id.vars = c('id', 'method', 'source', CONF_MAT_NAMES))

ggplot(per_melt) +
  geom_bar(
    mapping = aes(x = variable, y = value, fill = method),
    alpha = 0.8,
    stat = "identity",
    position = position_dodge2()
  ) +
  geom_text(
    aes(x = variable, y = -0.02, label = round(value, 2)),
    position = position_dodge2(width = 0.9),
    size = 3
  ) +
  facet_wrap(~source) +
  labs(title = "Performance Assessment", x = "metric", y = "performance") +
  scale_fill_brewer(palette = 'Set1') +
  #scale_fill_manual(name="Model", values=c("#66FF66", "#FF6600", "#CC0000", "#CC6699", "#9999CC"))+
  theme_classic() +
  theme(
    legend.position = 'top',
    axis.text.x = element_text(angle=45, vjust=1, hjust=1)
  )

ggsave(barplot_path, width = 8, height=8)


get_curve_dt <- function(predictions, curve) {
  if (curve == 'ROC') {
    curve_func <- roc_curve
    auc_func <- roc_auc
  } else if (curve == 'PR') {
    curve_func <- pr_curve
    auc_func <- pr_auc
  }
  curve_data <- lapply(
    predictions,
    function(dt) {
      dt %>%
        curve_func(truth = label, score, event_level = 'second') %>%
        arrange(.threshold) %>%
        as.data.table() -> res_dt
      auc_val <- auc_func(dt, label, score, event_level = 'second')
      res_dt[, auc := auc_val$.estimate]
      res_dt[, method := dt$method[1]]
      res_dt[, source := dt$source[1]]
      res_dt
    }
  )
  rbindlist(curve_data)
}


get_auc_anno <- function(curve_dt) {
  anno_dt <- unique(curve_dt[, .(method, source, auc)])
  anno_dt[, label := paste(method, ':', round(auc, 2))]
  anno_dt[, x := 0.75]
  anno_dt[, y := 0.07 * .I]
  anno_dt
}


# ROC curves
print('ROC curves')
roc_dt <- get_curve_dt(predictions, 'ROC')
roc_anno <- get_auc_anno(roc_dt)

ggplot(roc_dt, aes(1 - specificity, sensitivity, col = method)) +
  geom_path() +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") +
  facet_wrap(~source) +
  coord_equal() +
  geom_text(
      data = roc_anno,
        mapping = aes(x, y, col = method),
        label = roc_anno$label,
        inherit.aes = FALSE
  ) +
  labs(title = 'ROC Curves') +
  scale_color_brewer(palette = 'Dark2') +
  theme_classic() +
  theme(legend.position = 'top')

ggsave(roc_path, width = 8, height = 9)


# Precision-Recall curves
print('Precision-recall curve')
pr_dt <- get_curve_dt(predictions, 'PR')
pr_anno <- get_auc_anno(pr_dt)

ggplot(pr_dt, aes(recall, precision, col = method)) +
  geom_path() +
  facet_wrap(~source) +
  coord_equal() +
  geom_text(
    data = pr_anno,
    mapping = aes(x, y, col = method),
    label = pr_anno$label,
    inherit.aes = FALSE
  ) +
  labs(title = 'Precision-recall Curves') +
  scale_color_brewer(palette = 'Dark2') +
  theme_classic() +
  theme(legend.position = 'top')

ggsave(pr_path, width = 8, height = 9)

