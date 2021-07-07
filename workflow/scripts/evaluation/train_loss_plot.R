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

train_df <- fread(snakemake@input$train_eval)
train_long <- melt(train_df, id.vars = 'epoch', variable.name = 'metric')

ggplot(train_long, aes(epoch, value, col = metric)) +
  geom_line() +
  scale_color_brewer(palette = 'Set1') +
  theme_classic()

ggsave(snakemake@output$plot)
