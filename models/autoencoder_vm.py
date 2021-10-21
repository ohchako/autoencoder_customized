# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 15:38:21 2021

@author: Ouchi
"""
from tensorflow.keras.layers import Input, Dense, concatenate
from tensorflow.keras.models import Model
from keras import backend as K




def model():
    encoding_dim = 100
    input_dim = 2001
    output_dim = 2001


    Input1 = Input(shape=(input_dim,))
    Input2 = Input(shape=(input_dim,))
    Input3 = Input(shape=(input_dim,))
    Input4 = Input(shape=(input_dim,))
    Input5 = Input(shape=(input_dim,))

    Dense1 = Dense(encoding_dim,activation='relu')(Input1)
    Dense2 = Dense(encoding_dim,activation='relu')(Input2)
    Dense3 = Dense(encoding_dim,activation='relu')(Input3)
    Dense4 = Dense(encoding_dim,activation='relu')(Input4)
    Dense5 = Dense(encoding_dim,activation='relu')(Input5)

    Merge = concatenate([Dense1,Dense2,Dense3,Dense4,Dense5])

    ZZ = Dense(encoding_dim, activation='relu')(Merge)

# encodeされたデータを再構成した波形を格納する変数
    Decoded = Dense(output_dim,activation='sigmoid')(ZZ)

# 入力波形を再構成するModelとして定義
    model = Model(inputs=[Input1,Input2,Input3,Input4,Input5], outputs=[Decoded])
    
    # model.summary()

    # モデル可視化
    # from keras.utils import plot_model
    # plot_model(model, show_shapes = True)

    def root_mean_squared_error(y_true, y_pred):
            return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1)) 

    model.compile(optimizer='adam', loss = root_mean_squared_error, metrics = root_mean_squared_error)
    
    return model