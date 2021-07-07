import keras
import pandas as pd

if __name__ == "__main__":
    import sys

    sys.path.append(snakemake.input['common'])
    from common import convert_seqs_to_matrix

    model = keras.models.load_model(filepath=snakemake.input['model'])
    test_set = pd.read_table(snakemake.input['test_data'].__str__())
    x_test = convert_seqs_to_matrix(test_set['encoded_seq'])
    y_test = keras.utils.to_categorical(test_set['label'], 2)

    predictions = model.predict(x_test)
    df = pd.DataFrame(y_test[:, 1], columns=['label'])
    df['score'] = predictions[:,0].astype(int)
    df['prediction'] = predictions[:,1]

    pd.to_csv(df, snakemake.output['prediction'], sep='\t')
