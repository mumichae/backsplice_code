#DeepCirCode

## 1 libraries
library(keras)
library(ROCR)
library(reticulate)

if (!requireNamespace("DeepCirCode")) {
  devtools::install_github("BioDataLearning/DeepCirCode")
}
library(DeepCirCode)

#specify the python environment you want to work with, which has all properties of Tenserflow and Keras
# use_python(, required = TRUE)

#check of keras package
keras_model_sequential()

## Getting the raw training and test sets

data(HumanTrain)  # datasets are lazy loaded as variable “HumanTrain” until being used
data(HumanTest)

# Check datasets:

dim(HumanTrain)
colnames(HumanTrain)

HumanTrain$label[1]
HumanTrain$RNA_input_seq[1]
HumanTrain$encoded_seq[1]



##Converting one-hot encoded sequence into 3D matrix as the input for DeepCirCode

# Prepare x_train and x_test as 3D matrix:

x_train <- convStringToMatrix(HumanTrain$encoded_seq)
x_test <- convStringToMatrix(HumanTest$encoded_seq)

# This process may take several minutes
# x_train[1,2,] will print out the one-hot encoded vector of the 2nd nucleotide of HumanTrain$RNA_input_seq[1]



# Prepare y_train and y_test, use Keras function:

y_train <- to_categorical(HumanTrain$label,2)
y_test <- to_categorical(HumanTest$label,2)


##Testing DeepCirCode with test sets

data(x_test_human)
data(y_test_human)

##Applying DeepCirCode on the test sets:

DeepCirCode_test(x_test_human, y_test_human, "human") # don't forget the "" for species argument

# The DeepCirCode model architecture, confusion matrix, loss, accuracy, ROC curve and ROC_AUC value will be printed.
