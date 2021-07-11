import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import svm
import warnings
import joblib

if __name__ == "__main__":

    freq_train = pd.read_table(snakemake.input['train_features'].__str__()
                               , header=None).iloc[:, 1:]
    y_train = pd.read_table(
        snakemake.input['train_labels'].__str__(),
        dtype='string', header=None)
    y_train = y_train.iloc[:, 0].str.replace('01', '1')
    y_train = y_train.str.replace('10', '0')

    model_path = snakemake.output['model'].__str__()

    # we need to extract important sequence features!
    # frequencies of 1/2/3-mer compositions, normalized by length of their intron/exon -> 336 k-mer compositional features!!!



    # Fit the SVM model
    warnings.filterwarnings('ignore')
    model = svm.SVC(C=1, kernel='rbf', gamma=1, probability=True)
    model.fit(freq_train, y_train)

    # svm_model$SV

    #plot(svm_model, xy_train_human, most important vector, second most important vector, slice = list(?))

    joblib.dump(model, model_path)
