import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import svm
import warnings
import joblib


if __name__ == "__main__":
    freq_test = pd.read_table(snakemake.input['test_features'].__str__()
                               , header=None).iloc[:, 1:]
    y_test= pd.read_table(
        snakemake.input['test_labels'].__str__(),
        dtype='string', header=None)
    y_test = y_test.iloc[:, 0].str.replace('01', '1')
    y_test = y_test.str.replace('10', '0')

    model_path = snakemake.input['model'].__str__()
    prediction_path = snakemake.output['prediction'].__str__()
    output = open(prediction_path, "w")

    model = joblib.load(model_path)

    # Make prediction
    prediction = model.predict(freq_test)
    score = model.predict_proba(freq_test)[:, 1]

    # Get results
    result = pd.DataFrame()
    result['label'] = y_test
    result['score'] = score.tolist()
    result['prediction'] = prediction.tolist()
    print(result.head())

    result.to_csv(prediction_path, sep='\t', index=False)

