# Title     : Wang et al. 2019 dataset
# Objective : Load data
# Created by: mumichae
# Created on: 6/6/21
library(DeepCirCode)
library(data.table)
library(stringr)

# count sequence features of the 1/2/3-mers in the 4 flanking sequences
# get all combinations 1/2/3-mers (3-mers copied from a codon table)
combinations <- c(
  "A", "C", "G", "U",
  "AA", "AC", "AG", "AU",
  "CA", "CC", "CG", "CU",
  "GA", "GC", "GG", "GU",
  "UA", "UC", "UG", "UU",
  "AUG", "UGG", "UAU", "UAC", "UUU", "UUC", "UGU", "UGC",
  "AAU", "AAC", "GAU", "GAC", "CAA", "CAG", "GAA", "GAG",
  "CAU", "CAC", "AAA", "AAG", "AUU", "AUC", "AUA", "GGU",
  "GGC", "GGA", "GGG", "GCU", "GCC", "GCA", "GCG", "GUU",
  "GUC", "GUA", "GUG", "ACU", "ACC", "ACA", "ACG", "CCU",
  "CCC", "CCA", "CCG", "CUU", "CUC", "CUA", "CUG", "UUA",
  "UUG", "UCU", "UCC", "UCA", "UCG", "AGU", "AGC", "CGU",
  "CGC", "CGA", "CGG", "AGA", "AGG", "UAA", "UAG", "UGA"
)


get_frequencies <- function(dataset) {

  output <- data.frame(matrix(nrow = nrow(dataset), ncol = 0))
  names <- list()

  for (seq_part in 0:3) {

    print(seq_part)
    part <- substring(dataset[, 2], seq_part * 50 + 1, seq_part * 50 + 50)

    for (c in combinations) {

      column <- c()

      for (entity in part) {
        indices <- gregexpr(paste("(?=", c, ")", sep = ""), entity, perl = TRUE)[[1]]
        if (indices[1] == -1) { #that means no match
          frequency = 0
        }
        else {
          frequency = length(indices) / str_length(entity)
        }
        column <- c(column, frequency)
      }

      names = c(names, paste(seq_part, c, sep = "_"))
      output <- cbind(output, column)
    }

  }
  colnames(output) <- names
  return(output)
}

data(HumanTrain)  # datasets are lazy loaded as variable "HumanTrain" until being used
data(HumanTest)
data(y_train_human)
data(y_test_human)

freq_train <- get_frequencies(HumanTrain)
freq_test <- get_frequencies(HumanTest)

saveRDS(freq_train, snakemake@output$train)
saveRDS(freq_train, snakemake@output$test)

# TODO: save all necessary output
fwrite(HumanTrain, snakemake@output$train, sep = '\t')
