# Adapted from https://github.com/BioDataLearning/DeepCirCode.git
# inst/get_motifs.py
import keras.utils
import pandas as pd
from keras.models import Sequential
from keras.layers import Conv1D, MaxPooling1D
from keras.layers import Dense, Dropout, Flatten
from keras import losses

if __name__ == "__main__":
    import sys

    sys.path.append(snakemake.input['common'])
    from common import convert_seqs_to_matrix

    train_set = pd.read_table(snakemake.input['train_data'].__str__())
    x_train = convert_seqs_to_matrix(train_set['encoded_seq'])
    y_train = keras.utils.to_categorical(train_set['label'], 2)

    model = Sequential()

    conv_layer1 = Conv1D(filters=128,
                         kernel_size=12,
                         strides=1,
                         padding='valid',
                         activation='relu',
                         input_shape=(200, 4))

    conv_layer2 = Conv1D(filters=128,
                         kernel_size=6,
                         strides=1,
                         padding='valid',
                         activation='relu')

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
        optimizer="RMSprop",
        metrics=['accuracy']
    )

    print('fit model')
    model.fit(
        x_train, y_train,
        batch_size=128,
        epochs=80,
        validation_split=0.1
    )

    model.save(snakemake.output['model'].__str__())
