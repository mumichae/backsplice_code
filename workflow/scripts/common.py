# Any Python script in the scripts folder will be able to import from this module.
import numpy as np


def convert_seqs_to_matrix(seqs):
    x_array = np.zeros(shape=(len(seqs), 200, 4))
    for i, seq in enumerate(seqs):
        j = 0
        for c in range(200):
            for r in range(4):
                x_array[i, c, r] = seq[j]
                j += 1
    return x_array
