#DeepCirCode

## 1 libraries
library(keras)
library(ROCR)
library(reticulate)

#if (!requireNamespace("DeepCirCode")) {
#  devtools::install_github("BioDataLearning/DeepCirCode")
#}

library(DeepCirCode)

#specify the python environment you want to work with, which has all properties of Tenserflow and Keras
# use_python(, required = TRUE)

#check of keras package
keras_model_sequential()

## Getting the raw training and test sets

# data(HumanTrain)  # datasets are lazy loaded as variable “HumanTrain” until being used
# data(HumanTest)
test_set <- read.table("//wsl$/Ubuntu-20.04/home/laura/ML4RG_project/backsplice_code/results/processed_data/features/DeepCirCode/DiLiddo2019/all_data.tsv",
                       header = TRUE, colClasses = 'character')            #(snakemake@input$test_set)
train_set <- read.table("//wsl$/Ubuntu-20.04/home/laura/ML4RG_project/backsplice_code/results/processed_data/features/DeepCirCode/Wang2019/all_data.tsv",
                        header = TRUE, colClasses = 'character') 

prediction_path <- "//wsl$/Ubuntu-20.04/home/laura/ML4RG_project/backsplice_code/results/evaluation/DeepCirCode/Wang2019/prediction.tsv" # snakemake@output$prediction

# Check datasets:

dim(test_set)
colnames(test_set)

test_set$label[1]
test_set$RNA_input_seq[1]
test_set$encoded_seq[1]



##Converting one-hot encoded sequence into 3D matrix as the input for DeepCirCode

# Prepare x_train and x_test as 3D matrix:

x_train <- convStringToMatrix(train_set$encoded_seq)
x_test <- convStringToMatrix(test_set$encoded_seq)

# This process may take several minutes
# x_train[1,2,] will print out the one-hot encoded vector of the 2nd nucleotide of HumanTrain$RNA_input_seq[1]



# Prepare y_train and y_test, use Keras function:

y_train <- to_categorical(train_set$label,2)
y_test <- to_categorical(test_set$label,2)


##Testing DeepCirCode with test sets

# data(x_test_human)
# data(y_test_human)

##Applying DeepCirCode on the test sets:

result <- DeepCirCode_test(x_test, y_test, "human") # don't forget the "" for species argument

# The DeepCirCode model architecture, confusion matrix, loss, accuracy, ROC curve and ROC_AUC value will be printed.


#Part 2:
# train own model based on DCCs architecture.

path <- system.file("getMotifs.py", package = "DeepCirCode")
source_python(path) 

# Get position probability matrix for each kernel: 

DeepCirCode_getMotifs(x_train, y_train, x_test, y_test) 

# This step may take tens of minutes since training epochs = 80
# When finished, the "DeepCirCode_Position_Probability_Matrix.txt" file will be saved in **your current working directory** 
