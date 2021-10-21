# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:20:04 2021

for 2 cells prediction

@author: Ouchi
"""


from tensorflow.keras.layers import Input, Dense, concatenate
from tensorflow.keras.models import Model
from keras import backend as K


def model_2():
    encoding_dim = 100
    input_dim = 2001
    output_dim = 2001


    Input1 = Input(shape=(input_dim,))
    Input2 = Input(shape=(input_dim,))
    
    Dense1 = Dense(encoding_dim,activation='relu')(Input1)
    Dense2 = Dense(encoding_dim,activation='relu')(Input2)
    
    Merge = concatenate([Dense1,Dense2])

    ZZ = Dense(encoding_dim, activation='relu')(Merge)

# encodeされたデータを再構成した波形を格納する変数
    Decoded = Dense(output_dim,activation='sigmoid')(ZZ)

# 入力波形を再構成するModelとして定義
    model_2 = Model(inputs=[Input1,Input2], outputs=[Decoded])
    
    # model.summary()

    # モデル可視化
    # from keras.utils import plot_model
    # plot_model(model, show_shapes = True)

    def root_mean_squared_error(y_true, y_pred):
            return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1)) 

    model_2.compile(optimizer='adam', loss = root_mean_squared_error, metrics = root_mean_squared_error)
    
    return model_2