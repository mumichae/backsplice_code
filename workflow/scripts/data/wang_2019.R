# Extract the features for SVM and RF from the data of Wang 2019 or 
# our own training and test data set


#library(DeepCirCode)
library(data.table)
library(stringr)

# count sequence features of the 1/2/3-mers in the 4 flanking sequences
# get all combinations 1/2/3-mers (3-mers copied from a codon table)

# regex of combinations
combinations <- c(
  "A", "C", "G", "[UT]",
  "AA", "AC", "AG", "A[UT]",
  "CA", "CC", "CG", "C[UT]",
  "GA", "GC", "GG", "G[UT]",
  "[UT]A", "[UT]C", "[UT]G", "[UT][UT]",
  "A[UT]G", "[UT]GG", "[UT]A[UT]", "[UT]AC", "[UT][UT][UT]", "[UT][UT]C", "[UT]G[UT]", "[UT]GC",
  "AA[UT]", "AAC", "GA[UT]", "GAC", "CAA", "CAG", "GAA", "GAG",
  "CA[UT]", "CAC", "AAA", "AAG", "A[UT][UT]", "A[UT]C", "A[UT]A", "GG[UT]",
  "GGC", "GGA", "GGG", "GC[UT]", "GCC", "GCA", "GCG", "G[UT][UT]",
  "G[UT]C", "G[UT]A", "G[UT]G", "AC[UT]", "ACC", "ACA", "ACG", "CC[UT]",
  "CCC", "CCA", "CCG", "C[UT][UT]", "C[UT]C", "C[UT]A", "C[UT]G", "[UT][UT]A",
  "[UT][UT]G", "[UT]C[UT]", "[UT]CC", "[UT]CA", "[UT]CG", "AG[UT]", "AGC", "CG[UT]",
  "CGC", "CGA", "CGG", "AGA", "AGG", "[UT]AA", "[UT]AG", "[UT]GA"
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
saveRDS(freq_test, snakemake@output$test)

# TODO: save all necessary output
#fwrite(HumanTrain, snakemake@output$train, sep = '\t')
