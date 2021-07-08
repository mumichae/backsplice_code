# Title     : Training Performance Plots
# Objective : Plot diagnostic plots
# Created by: mumichae
# Created on: 7/7/21

saveRDS(snakemake, ".snakemake/loss_plot.RDS")

# assess performance of our models
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

train_dt <- fread(snakemake@input$train_stats)
train_long <- melt(train_dt, id.vars = c('epoch', 'dataset'), variable.name = 'metric')

ggplot(train_long, aes(epoch, value, col = dataset)) +
  geom_line() +
  facet_wrap(~metric) +
  scale_color_brewer(palette = 'Set1') +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(snakemake@output$plot)
