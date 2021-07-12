import umap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    
    import sys
    sys.path.append(snakemake.input['common'])
    from common import convert_seqs_to_matrix
    
    train_data = pd.read_table(snakemake.input['train'])
    test_data = pd.read_table(snakemake.input['test'])
    
    train_data['data'] = 'train'
    test_data['data'] = 'test'
    
    all_data = pd.concat([train_data, test_data])
    
    print('Encode sequences')
    X = convert_seqs_to_matrix(all_data['encoded_seq'])
    X = X.reshape((len(X), -1))
    
    print('UMAP embedding')
    reducer = umap.UMAP(n_components=2)
    emb = reducer.fit_transform(X)
    
    figure, (ax1, ax2) = plt.subplots(1, 2)

    ax1.scatter(
        emb[:, 0],
        emb[:, 1],
        c=['navy' if x == 'train' else 'yellow' for x in all_data['data']],
        alpha=0.1,
        s=1
    )
    ax1.set_title('By test/train')

    ax2.scatter(
        emb[:, 0],
        emb[:, 1],
        c=all_data['label'],
        alpha=0.1,
        s=1
    )
    ax2.set_title('By labels')

    plt.subplots_adjust(hspace = 0.5)
    figure.savefig(snakemake.output['plot'], dpi=200)
