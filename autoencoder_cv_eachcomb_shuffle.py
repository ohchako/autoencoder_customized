# -*- coding: utf-8 -*-
"""
Created on Tue Jun 29 11:09:23 2021

perform cross-validation --shuffle(combination)--

@author: Ouchi
"""


import numpy as np

import os

import itertools

from models.autoencoder_vm1 import model_1
from models.autoencoder_vm2 import model_2
from models.autoencoder_vm3 import model_3
#from models.autoencoder_vm4 import model_4


path = 'E:\\vitro_autoencoder_data\\210701\\3cells\\190910'
os.chdir(path)
# os.makedirs('shuffle_eachcomb_cv')


# cellの数【データの記録細胞数が異なれば変更】
cellnum = 3
# dataをimportする(波形ごとに[0 1]でスケーリング前処理済み)→ x
data = np.loadtxt('190910_3cells_normalized.txt',delimiter =',',dtype = float)
# 膜電位dataをimportする(波形ごとに[0 1]でスケーリング前処理済み)→ y
Vm = np.loadtxt('190910_vmtrace_normalized_3cells.txt',delimiter =',',dtype = None)

index_list = np.array(range(data.shape[0]))

# matlabでrandperm関数で作成（インデックスが1からなので注意）→210319に乱数の取り込み方法をアップデートした
P = np.loadtxt('190910_Vmshuffleno.txt',delimiter =',',dtype = 'int32')

def split_list(l,i):
    for idx in range(0,l.shape[1],i):
        yield l[:,idx:idx+i]


Vm2 = np.array(list(split_list(Vm,2001)))

split_random = list(split_list(P,100))
p1 = split_random[0]
p2 = split_random[1]
p3 = split_random[2]
# p4 = split_random[3]


###--- cross-validation
# 5回のcvですべてのSWが一回はテストデータとなるようにする
cv_times = 5
# ★ここはデータごとにしばらく要確認！！！！！！！！！！　２１０７０１21:43
block_s = data.shape[0] // cv_times
if (block_s * cv_times) < data.shape[0]:
    block_s = block_s + 1

## 各組合せごとに10回シャッフル           
for z in range(1,cellnum):
    if z == 2:
        break#組合せごとのファイルを作成
    
    a = np.arange(cellnum)
    A = list(itertools.combinations(a,z))
    dirname = str(cellnum) +'C' + str(z)
    path_s = path +'\\shuffle_eachcomb_cv'
    os.chdir(path_s)
    os.makedirs(dirname)
    os.chdir(dirname)
    os.makedirs('x_test')
    os.makedirs('decoded_vm')
    path_x_test = path_s + '\\' + dirname + '\\x_test'
    path_decoded_vm = path_s + '\\' + dirname + '\\decoded_vm'
    
    for y in range(len(A)):
        p = A[y]
        Vm3 = Vm2[p,:,:]#組合せごとに必要なcellを取り出す
        Vm4 = np.hstack(Vm3) #扱いやすいので一時的に２次元にする
    
        for x in range(100):
           
            n = x
            if Vm3.shape[0] == 1:
                P1 = p1[:,n]-1
                a2 = Vm3[0,P1]
                Vm5 = a2
            elif Vm3.shape[0] == 2:
                P1 = p1[:,n]-1
                a2 = Vm3[0,P1]
                P2 = p2[:,n]-1
                b2 = Vm3[1,P2]
                Vm5 = np.concatenate([a2,b2],1)
            elif Vm3.shape[0] == 3:
                P1 = p1[:,n]-1
                a2 = Vm3[0,P1]
                P2 = p2[:,n]-1
                b2 = Vm3[1,P2]
                P3 = p3[:,n]-1
                c2 = Vm3[2,P3]
                Vm5 = np.concatenate([a2,b2,c2],1)
            # elif Vm3.shape[0] == 4:
            #     P1 = p1[:,n]-1
            #     a2 = Vm3[0,P1]
            #     P2 = p2[:,n]-1
            #     b2 = Vm3[1,P2]
            #     P3 = p3[:,n]-1
            #     c2 = Vm3[2,P3]
            #     P4 = p4[:,n]-1
            #     d2 = Vm3[3,P3]
            #     Vm5 = np.concatenate([a2,b2,c2,d2],1)
            
            for cv in range(cv_times):
                idx = np.arange(block_s) + ((block_s) * cv)
                if cv == int(cv_times - 1):
                    idx = np.arange(block_s) + ((block_s) * cv) - (((block_s) * (cv + 1)) - data.shape[0])
            
                x_test = data[idx]
                y_test = Vm5[idx]
                y_test_vms = np.array(list(split_list(y_test,2001)))
        
        # testに割り当てた以外のidxの9割をtraining用にする
                idx_d = np.delete(index_list,idx)
                x_train = data[idx_d]
                y_train = Vm5[idx_d]
                y_train_vms = np.array(list(split_list(y_train,2001)))
            
                if Vm3.shape[0] == 1:
                    autoencoder = model_1()
        
                    autoencoder.fit([y_train_vms[0,:,:]], x_train,
                                    epochs=100,
                                    batch_size=16,
                                    shuffle=True,
                                    validation_data=([y_test_vms[0,:,:]], x_test))
                
                    predicted = autoencoder.predict({"Input1":y_test_vms[0,:,:]})
            
                elif Vm3.shape[0] == 2:
                    autoencoder = model_2()
                
                    autoencoder.fit([y_train_vms[0,:,:],y_train_vms[1,:,:]], x_train,
                                    epochs=100,
                                    batch_size=16,
                                    shuffle=True,
                                    validation_data=([y_test_vms[0,:,:],y_test_vms[1,:,:]], x_test))
                    predicted = autoencoder.predict({"Input1":y_test_vms[0,:,:],"Input2":y_test_vms[1,:,:]})
            
                elif Vm3.shape[0] == 3:
                    autoencoder = model_3()
            
                    autoencoder.fit([y_train_vms[0,:,:],y_train_vms[1,:,:],y_train_vms[2,:,:]], x_train,
                                    epochs=100,
                                    batch_size=16,
                                    shuffle=True,
                                    validation_data=([y_test_vms[0,:,:],y_test_vms[1,:,:],y_test_vms[2,:,:]], x_test))
                    predicted = autoencoder.predict({"Input1":y_test_vms[0,:,:],"Input2":y_test_vms[1,:,:],"Input3":y_test_vms[2,:,:]})
            
                # elif Vm3.shape[0] == 4:
                #     autoencoder = model_4()
                
                #     autoencoder.fit([y_train_vms[0,:,:],y_train_vms[1,:,:],y_train_vms[2,:,:],y_train_vms[3,:,:]], x_train,
                #                     epochs=100,
                #                     batch_size=16,
                #                     shuffle=True,
                #                     validation_data=([y_test_vms[0,:,:],y_test_vms[1,:,:],y_test_vms[2,:,:],y_test_vms[3,:,:]], x_test))
                #     predicted = autoencoder.predict({"Input1":y_test_vms[0,:,:],"Input2":y_test_vms[1,:,:],"Input3":y_test_vms[2,:,:],"Input4":y_test_vms[3,:,:]})
            
                # autoencoder.save(os.path.join(savepath, "model_{}.h5".format(cv)))
                # autoencoder.save(f"model_{z,y,x,cv}.h5".format(cv))
                os.chdir(path_decoded_vm)
                np.savetxt(f"combination_cv-{z,y,x,cv}.csv",predicted,fmt='%.6f',delimiter=',')
                os.chdir(path_x_test)
                np.savetxt(f"combination_cv_x_test-{z,y,x,cv}.csv",x_test,fmt='%.6f',delimiter=',')
                del autoencoder,y_test,y_train,x_test,idx,idx_d,y_test_vms,y_train_vms

