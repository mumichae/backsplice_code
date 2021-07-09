# Adapted from https://github.com/BioDataLearning/DeepCirCode.git
# inst/get_motifs.py
from tensorflow.python import keras
import keras.utils
import pandas as pd
from keras.models import Sequential
from keras.layers import Conv1D, MaxPooling1D
from keras.layers import Dense, Dropout, Flatten
from keras.regularizers import l2
from keras import losses
import numpy as np
from sklearn.model_selection import train_test_split

if __name__ == "__main__":
    import sys

    sys.path.append(snakemake.input['common'])
    from common import convert_seqs_to_matrix

    np.random.seed(123)

    train_set = pd.read_table(snakemake.input['train_data'].__str__())
    x_train = convert_seqs_to_matrix(train_set['encoded_seq'])
    y_train = keras.utils.to_categorical(train_set['label'], 2)
    x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.2)

    print(f'train: {len(y_train)}\t validation: {len(y_val)}')

    model = Sequential()

    conv_layer1 = Conv1D(filters=128,
                         kernel_size=12,
                         strides=1,
                         padding='valid',
                         activation='relu',
                         input_shape=(200, 4),
                         kernel_regularizer=l2(0.01)
                )

    conv_layer2 = Conv1D(filters=128,
                         kernel_size=6,
                         strides=1,
                         padding='valid',
                         activation='relu',
                         kernel_regularizer=l2(0.01)
                 )

    model.add(conv_layer1)
    model.add(Dropout(0))

    model.add(conv_layer2)
    model.add(Dropout(0.5))

    model.add(MaxPooling1D(pool_size=4, strides=4))
    model.add(Dropout(0.5))

    model.add(Flatten())
    model.add(Dense(2, activation='softmax'))

    model.summary()

    model.compile(
        loss=losses.binary_crossentropy,
        optimizer=keras.optimizers.RMSprop(learning_rate=1e-3),
        metrics=['accuracy', keras.metrics.Precision(), keras.metrics.Recall()]
    )

    print('fit model')
    history = model.fit(
        x_train, y_train,
        batch_size=128,
        epochs=80,
        validation_data=(x_val, y_val),
        callbacks=[
            keras.callbacks.EarlyStopping(
                monitor="val_loss",
                min_delta=1e-3,
                patience=3,
                verbose=1,
            )
        ]
    )

    df = pd.DataFrame.from_dict(history.history)
    n_epochs = df.shape[0]

    # long format
    train_cols = [x for x in df.columns if not x.startswith("val")]
    val_df = df[[x for x in df.columns if x.startswith("val")]]
    val_df.columns = train_cols
    df = pd.concat([df[train_cols], val_df])

    df['dataset'] = [*n_epochs*['train'], *n_epochs*['validation']]
    # add epochs
    df['epoch'] = 2*[*range(n_epochs)]

    df.to_csv(
        snakemake.output['training_stats'].__str__(),
        index=False,
        sep='\t'
    )
    model.save(snakemake.output['model'].__str__())
