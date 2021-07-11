# Adapted from https://github.com/BioDataLearning/DeepCirCode.git
# inst/get_motifs.py

from tensorflow.python import keras
import keras.utils
from keras import backend as K
import pandas as pd
import numpy as np


def get_motifs(model, prob_matrix_file):
    # motif visualization:
    ### motif visualization for BS model human:
    conv_output = model.get_layer('conv1d').get_output_at(0)
    f = K.function([model.input], [K.argmax(conv_output, axis=1), K.max(conv_output, axis=1)])

    n_filters = 128//4
    motifs = np.zeros((n_filters, 12, 4))
    nsites = np.zeros(n_filters)

    # select the positive samples from test set:
    row = y_test[:, 1] == 1.0
    x_test_positive = x_test[np.ix_(row)]

    for i in range(0, len(x_test_positive), 100):
        x = x_test_positive[i:i + 100]
        z = f([x])
        max_inds = z[0]  # N x M matrix, where M is the number of motifs
        max_acts = z[1]
        for m in range(n_filters):
            for n in range(len(x)):
                # Forward strand
                if max_acts[n, m] > 0:
                    nsites[m] += 1
                    motifs[m] += x[n, max_inds[n, m]:max_inds[n, m] + 12, :]

    print('Making motifs')
    motifs = motifs[:, :, [0, 3, 2, 1]]
    motifs_file = open(prob_matrix_file, 'w')
    motifs_file.write(
        'MEME version 4.9.0\n\n'
        'ALPHABET= ACGU\n\n'
        'strands: + -\n\n'
        'Background letter frequencies (from uniform background):\n'
        'A 0.25000 C 0.25000 G 0.25000 U 0.25000\n\n'
    )

    for m in range(n_filters):
        if nsites[m] == 0:
            continue
        motifs_file.write('MOTIF M_n%i O%i\n' % (m, m))
        motifs_file.write("letter-probability matrix: alength= 4 w= %i nsites= %i E= 1337.0e-6\n" % (12, nsites[m]))
        for j in range(12):
            motifs_file.write("%f %f %f %f\n" % tuple(1.0 * motifs[m, j, 0:4] / np.sum(motifs[m, j, 0:4])))
        motifs_file.write('\n')

    motifs_file.close()
    print("Done")


if __name__ == "__main__":
    import sys

    sys.path.append(snakemake.input['common'])
    from common import convert_seqs_to_matrix

    model = keras.models.load_model(filepath=snakemake.input['model'])
    test_set = pd.read_table(snakemake.input['test_data'].__str__())
    x_test = convert_seqs_to_matrix(test_set['encoded_seq'])
    y_test = keras.utils.to_categorical(test_set['label'], 2)

    preds = model.predict(x_test)
    res = model.evaluate(x_test, y_test)
    print(res)
    pd.DataFrame.from_dict(
        {
            'label': test_set['label'],
            'score': preds[:,1],
            'prediction': preds.argmax(axis=1)
        }
    ).to_csv(snakemake.output['prediction'], sep='\t', index=False)
    
    model.summary()
    get_motifs(model, snakemake.output['motifs'])
