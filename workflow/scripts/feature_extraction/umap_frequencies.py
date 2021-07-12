import umap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    train_data = pd.read_table(snakemake.input['train'].__str__(), header=None, index_col=0)
    test_data = pd.read_table(snakemake.input['test'].__str__(), header=None, index_col=0)
    train_labels = [int(l.strip()[0]) for l in open(snakemake.input['train_labels'].__str__(), 'r')]
    test_labels = [int(l.strip()[0]) for l in open(snakemake.input['test_labels'].__str__(), 'r')]
    
    train_data['data'] = 'train'
    test_data['data'] = 'test'
    train_data['label'] = train_labels
    test_data['label'] = test_labels
    
    all_data = pd.concat([train_data, test_data])
    X = all_data.iloc[:,:-2].to_numpy()
    
    print('UMAP embedding')
    reducer = umap.UMAP(n_components=2)
    emb = reducer.fit_transform(X[:,:-2])
    
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

    figure.savefig(snakemake.output['plot'], dpi=200)
