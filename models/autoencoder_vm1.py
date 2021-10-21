# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:21:06 2021

for 1 cell prediction

@author: Ouchi
"""

from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.models import Model
from keras import backend as K


def model_1():
    encoding_dim = 100
    input_dim = 4001
    output_dim = 4001


    Input1 = Input(shape=(input_dim,))
    
    Dense1 = Dense(encoding_dim,activation='relu')(Input1)


# encodeされたデータを再構成した波形を格納する変数
    Decoded = Dense(output_dim,activation='sigmoid')(Dense1)

# 入力波形を再構成するModelとして定義
    model_1 = Model(inputs=[Input1], outputs=[Decoded])
    
    # model.summary()

    # モデル可視化
    # from keras.utils import plot_model
    # plot_model(model, show_shapes = True)

    def root_mean_squared_error(y_true, y_pred):
            return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1)) 

    model_1.compile(optimizer='adam', loss = root_mean_squared_error, metrics = root_mean_squared_error)
    
    return model_1