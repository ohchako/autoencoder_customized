# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 18:53:53 2021

@author: Ouchi

updated; 210629 
5細胞記録なら5細胞の完全シャッフル
perform cross validation --shuffle--


"""

import numpy as np
import os

# from models.autoencoder_vm import model # for 5cells
# from models.autoencoder_vm4 import model_4 # for 4cells
from models.autoencoder_vm3 import model_3 # for 3cells
 
def split_list(l,i):
    for idx in range(0,l.shape[1],i):
        yield l[:,idx:idx+i]
        
path = 'E:\\vitro_autoencoder_data\\210701\\3cells\\190927_3'
os.chdir(path)
os.makedirs('shuffle_cv')

# dataをimportする(波形ごとに[0 1]でスケーリング前処理済み)→ x
data = np.loadtxt('190927_3_3cells_normalized.txt',delimiter =',',dtype = float)
# 膜電位dataをimportする(波形ごとに[0 1]でスケーリング前処理済み)→ y
Vm = np.loadtxt('190927_3_vmtrace_normalized_3cells.txt',delimiter =',',dtype = None)

index_list = np.array(range(data.shape[0]))

split_vm = np.array(list(split_list(Vm,2001)))

# matlabでrandperm関数で作成（インデックスが1からなので注意）
P = np.loadtxt('190927_3_Vmshuffleno.txt',delimiter =',',dtype = 'int32')

split_random = list(split_list(P,100))
p1 = split_random[0]
p2 = split_random[1]
p3 = split_random[2]
# p4 = split_random[3]
# p5 = split_random[4]

###--- cross-validation
# 5回のcvですべてのSWが一回はテストデータとなるようにする
cv_times = 5
block_s = data.shape[0] // cv_times
if (block_s * cv_times) < data.shape[0]:
    block_s = block_s + 1

for x in range(10):
    path_s = path +'\\shuffle_cv'
    os.chdir(path_s)
    n = x
    P1 = p1[:,n]-1
    a2 = split_vm[0,P1]
    P2 = p2[:,n]-1
    b2 = split_vm[1,P2]
    P3 = p3[:,n]-1
    c2 = split_vm[2,P3]
    # P4 = p4[:,n]-1
    # d2 = split_vm[3,P4]
    # P5 = p5[:,n]-1
    # e2 = split_vm[4,P5]
    Vm_shuffle = np.concatenate([a2,b2,c2],1) # 乱数表に従って順番がシャッフルされた膜電位
    
    for cv in range(cv_times):
        idx = np.arange(block_s) + ((block_s) * cv)
        if cv == int(cv_times - 1):
            idx = np.arange(block_s) + ((block_s) * cv) - (((block_s) * (cv + 1)) - data.shape[0])
            
        x_test = data[idx]
        y_test = Vm_shuffle[idx]
        y_test_split = list(split_list(y_test,2001))
        y_test_vms = np.array(y_test_split)
        
        # testに割り当てた以外のidxの9割をtraining用にする
        idx_d = np.delete(index_list,idx)
        x_train = data[idx_d]
        y_train = Vm_shuffle[idx_d]
        y_train_split = list(split_list(y_train,2001))
        y_train_vms = np.array(y_train_split)
        
        autoencoder = model_3()

        autoencoder.fit(
        [y_train_vms[0,:,:],y_train_vms[1,:,:],y_train_vms[2,:,:]],
        x_train,
        epochs=100,
        batch_size=16,
        shuffle=True,
        validation_data=([y_test_vms[0,:,:],y_test_vms[1,:,:],y_test_vms[2,:,:]], x_test),
        )
        
        # autoencoder.save(os.path.join(path_s, f"model_{cv,x}.h5".format(cv)))
        
        predicted = autoencoder.predict({"Input1":y_test_vms[0,:,:],"Input2":y_test_vms[1,:,:],"Input3":y_test_vms[2,:,:]})
        np.savetxt(f"shuffle-{cv,x}.csv",predicted,fmt='%.6f',delimiter=',')
        np.savetxt(f"x_test-{cv,x}.csv",x_test,fmt='%.6f',delimiter=',')
        del x_test,y_test,y_train,y_test_split,y_test_vms,idx,idx_d,y_train_split,y_train_vms
    del n


# from sklearn.metrics import mean_squared_error # モデル評価用(平均二乗誤差)

# rmse = []
# rmse_all = []
# for x in range(len(x_test)):
#     mse = mean_squared_error(x_test[x,:], Decoded_img[x,:]) # MSE(平均二乗誤差)の算出
#     rmse = np.sqrt(mse) # RSME = √MSEの算出
#     rmse_all = np.hstack([rmse_all,rmse])

# rmse_mean = np.mean(rmse_all)
# print('RMSE :',rmse_mean)
