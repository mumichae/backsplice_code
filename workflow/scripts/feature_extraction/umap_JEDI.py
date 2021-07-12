import umap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def pad(x):
    y = np.zeros(4 * 128, dtype=np.int32)
    RL = min(len(x), 4 * 128)
    y[:RL] = x[:RL]
    return y


def flatten_features(series, data, label):
    df = pd.DataFrame.from_records(series.apply(lambda x: np.ravel(x)).apply(pad).values)
    df['data'] = data
    df['label'] = label
    return df


if __name__ == "__main__":

    train_data = pd.read_json(snakemake.input.train.__str__(), lines=True)
    test_data = pd.read_json(snakemake.input.test.__str__(), lines=True)
    
    all_acceptors = pd.concat(
        [
            flatten_features(train_data['acceptors'], 'train', train_data['label']),
            flatten_features(test_data['acceptors'], 'test', test_data['label'])
        ]
    )
    X_acc = all_acceptors.values[:, :-2]

    all_donors = pd.concat(
        [
            flatten_features(train_data['donors'], 'train', train_data['label']),
            flatten_features(test_data['donors'], 'test', test_data['label'])
        ]
    )
    X_don = all_donors.values[:, :-2]
    
    reducer = umap.UMAP(n_components=2)
    
    print('embed acceptors')
    emb_acc = reducer.fit_transform(X_acc)
    print('embed donors')
    emb_don = reducer.fit_transform(X_don)
    
    figure, axis = plt.subplots(2, 2)

    axis[0, 0].scatter(
        emb_acc[:, 0],
        emb_acc[:, 1],
        c=['navy' if x == 'train' else 'yellow' for x in all_acceptors['data']],
        alpha=0.1,
        s=1
    )
    axis[0, 0].set_title('Acceptors by test/train')

    axis[0, 1].scatter(
        emb_acc[:, 0],
        emb_acc[:, 1],
        c=all_acceptors['label'],
        alpha=0.1,
        s=1
    )
    axis[0, 1].set_title('Acceptors by labels')

    axis[1, 0].scatter(
        emb_don[:, 0],
        emb_don[:, 1],
        c=['navy' if x == 'train' else 'yellow' for x in all_donors['data']],
        alpha=0.1,
        s=1
    )
    axis[1, 0].set_title('Donors by test/train')

    axis[1, 1].scatter(
        emb_don[:, 0],
        emb_don[:, 1],
        c=all_donors['label'],
        alpha=0.1,
        s=1
    )
    axis[1, 1].set_title('Donors by labels')
    
    plt.subplots_adjust(hspace = 0.5)
    figure.savefig(snakemake.output['plot'], dpi=200)
