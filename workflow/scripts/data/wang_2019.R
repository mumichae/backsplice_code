# Extract the features for SVM and RF from the data of Wang 2019 or 
# our own training and test data set


#library(DeepCirCode)
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

#train <- data(HumanTrain)  # datasets are lazy loaded as variable "HumanTrain" until being used
#test <- data(HumanTest)
#data(y_train_human)
#data(y_test_human)

train <- data.frame(read.table(snakemake@input$train_data, sep = "\t", header = TRUE))
test <-  data.frame(read.table(snakemake@input$test_data, sep = "\t", header = TRUE))

#train <- data.frame(read.table('//wsl$/Ubuntu-20.04/home/laura/ML4RG_project/backsplice_code/results/processed_data/features/DeepCirCode/DiLiddo2019/all_data.tsv', sep = "\t", header = TRUE))


freq_train <- get_frequencies(train)
freq_test <- get_frequencies(test)

saveRDS(freq_train, snakemake@output$train)
saveRDS(freq_train, snakemake@output$test)

# TODO: save all necessary output
#fwrite(HumanTrain, snakemake@output$train, sep = '\t')
