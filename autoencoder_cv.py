# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 16:16:39 2021

perform cross-validation

@author: Ouchi
"""


import numpy as np
import os
#import matplotlib.pyplot as plt
#from pylab import rcParams

from models.autoencoder_vm import model
from models.autoencoder_vm4 import model_4 # for 4cells
from models.autoencoder_vm3 import model_3 # for 3cells
from models.autoencoder_vm1 import model_1 # for 1cell


savepath = 'F:\\vivo patch解析\\MC\\180219'
os.chdir(savepath)
os.makedirs('real_cv')


##---- 
# dataをimportする(波形ごとに[0 1]でスケーリング前処理済み)→ x
data = np.loadtxt('180219_LFP4001_n.txt',delimiter =',',dtype = float)
# 膜電位dataをimportする(波形ごとに[0 1]でスケーリング前処理済み)→ y（5cellsならSW数x10005）
Vm = np.loadtxt('180219_vm4001_n.txt',delimiter =',',dtype = None)

index_list = np.array(range(data.shape[0]))

def split_list(l,i):
    for idx in range(0,l.shape[1],i):
        yield l[:,idx:idx+i]
        

os.chdir('real_cv')

###--- cross-validation
# 5回のcvですべてのSWが一回はテストデータとなるようにする
cv_times = 5
block_s = data.shape[0] // cv_times
if (block_s * cv_times) < data.shape[0]:
    block_s = block_s + 1

for cv in range(cv_times):
    idx = np.arange(block_s) + ((block_s) * cv)
    if cv == int(cv_times - 1):
        idx = np.arange(block_s) + ((block_s) * cv) - (((block_s) * (cv + 1)) - data.shape[0])

    x_test = data[idx]
    y_test = Vm[idx]
    y_test_split = list(split_list(y_test,4001))
    y_test_vms = np.array(y_test_split)
        
        # testに割り当てた以外のidxの9割をtraining用にする
    idx_d = np.delete(index_list,idx)
    x_train = data[idx_d]
    y_train = Vm[idx_d]
    y_train_split = list(split_list(y_train,4001))
    y_train_vms = np.array(y_train_split)
        
    if y_test_vms.shape[0] == 3:
        autoencoder = model_3()
        
        autoencoder.fit(
        [y_train_vms[0,:,:],y_train_vms[1,:,:],y_train_vms[2,:,:]],
        x_train,
        epochs=100,
        batch_size=16,
        shuffle=True,
        validation_data=([y_test_vms[0,:,:],y_test_vms[1,:,:],y_test_vms[2,:,:]], x_test),
        )
        predicted = autoencoder.predict({"Input1":y_test_vms[0,:,:],"Input2":y_test_vms[1,:,:],"Input3":y_test_vms[2,:,:]})
    
    elif y_test_vms.shape[0] == 4:
        autoencoder = model_4()
        
        autoencoder.fit(
        [y_train_vms[0,:,:],y_train_vms[1,:,:],y_train_vms[2,:,:],y_train_vms[3,:,:]],
        x_train,
        epochs=100,
        batch_size=16,
        shuffle=True,
        validation_data=([y_test_vms[0,:,:],y_test_vms[1,:,:],y_test_vms[2,:,:],y_test_vms[3,:,:]], x_test),
        )
        predicted = autoencoder.predict({"Input1":y_test_vms[0,:,:],"Input2":y_test_vms[1,:,:],"Input3":y_test_vms[2,:,:],"Input4":y_test_vms[3,:,:]})
        
    elif y_test_vms.shape[0] == 5:
        autoencoder = model()
        
        autoencoder.fit(
        [y_train_vms[0,:,:],y_train_vms[1,:,:],y_train_vms[2,:,:],y_train_vms[3,:,:],y_train_vms[4,:,:]],
        x_train,
        epochs=100,
        batch_size=16,
        shuffle=True,
        validation_data=([y_test_vms[0,:,:],y_test_vms[1,:,:],y_test_vms[2,:,:],y_test_vms[3,:,:],y_test_vms[4,:,:]], x_test),
        )
        predicted = autoencoder.predict({"Input1":y_test_vms[0,:,:],"Input2":y_test_vms[1,:,:],"Input3":y_test_vms[2,:,:],"Input4":y_test_vms[3,:,:],"Input5":y_test_vms[4,:,:]})
        
    elif y_test_vms.shape[0] == 1:
        autoencoder = model_1()
        
        autoencoder.fit(
        [y_train_vms[0,:,:]],
        x_train,
        epochs=100,
        batch_size=16,
        shuffle=True,
        validation_data=([y_test_vms[0,:,:]], x_test),
        )
        predicted = autoencoder.predict({"Input1":y_test_vms[0,:,:]})
            
        
            
    autoencoder.save(os.path.join(savepath, "model_{}.h5".format(cv)))

    np.savetxt(f"real-{cv}.csv",predicted,fmt='%.6f',delimiter=',')
    np.savetxt(f"x_test-{cv}.csv",x_test,fmt='%.6f',delimiter=',')
    del x_test,y_test,y_test_split,y_test_vms,idx_d,x_train,y_train,y_train_split,y_train_vms,predicted


## load model and reconstruct data
from keras.models import load_model
from keras import backend as K
	
def root_mean_squared_error(y_true, y_pred):
            return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1)) 
model = load_model('model_4.h5',custom_objects={'root_mean_squared_error': root_mean_squared_error})		


# N = 6
# # figの大きさいじれる
# rcParams['figure.figsize'] = 10,10
# plt.figure()
# plt.subplot(3,1,1)
# plt.plot(x_test[N,:], color = "black")
# plt.subplot(3,1,2)
# plt.plot(Decoded_img[N,:], color = "red") 
# plt.subplot(3,1,3)
# plt.plot(x_test[N,:], color = "black")
# plt.plot(Decoded_img[N,:], color = "red")

# from sklearn.metrics import mean_squared_error # モデル評価用(平均二乗誤差)

# rmse = []
# rmse_all = []
# for x in range(len(x_test)):
#     mse = mean_squared_error(x_test[x,:], Decoded_img[x,:]) # MSE(平均二乗誤差)の算出
#     rmse = np.sqrt(mse) # RSME = √MSEの算出
#     rmse_all = np.hstack([rmse_all,rmse])

# rmse_mean = np.mean(rmse_all)
# print('RMSE :',rmse_mean)
