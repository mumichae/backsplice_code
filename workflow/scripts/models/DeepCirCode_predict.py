from tensorflow.python import keras
import keras.utils
import pandas as pd
import numpy as np

if __name__ == "__main__":
    import sys

    sys.path.append(snakemake.input['common'])
    from common import convert_seqs_to_matrix

    model = keras.models.load_model(filepath=snakemake.input['model'])
    test_set = pd.read_table(snakemake.input['test_data'].__str__())
    x_test = convert_seqs_to_matrix(test_set['encoded_seq'])
    y_test = keras.utils.to_categorical(test_set['label'], 2)

    preds = model.predict(x_test)
    pd.DataFrame.from_dict(
        {
            'label': y_test[:,1],
            'score': preds[:,1],
            'prediction': np.round(preds[:,1]).astype(int)
        }
    ).to_csv(snakemake.output['prediction'], sep='\t', index=False)

