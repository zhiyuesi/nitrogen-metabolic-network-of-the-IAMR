# -*- coding: utf-8 -*-
'''
@Author  :  月司
@Email   :  815603884@qq.com
@File    :  Plot_the_error_of_the_linear_solution_model.py
@version :  1.0
@Desc    : 给标准的莫诺特方程来添加高斯噪声，然后使用线性模型进行模拟，绘制了在不同误差下的动图。
'''
import numpy as np
from math import log
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
# 刻度设定器
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FuncFormatter,IndexLocator,LinearLocator,AutoLocator

# 自己写的一些matplotlib的配置，
from matplotlib_config import *
# 进行线性拟合
import numpy as np
from sklearn import linear_model
from collections import defaultdict


# 读入SOUR的动力学数据
result_path="./data"
# 模拟的高斯噪声数据分析
no_noise_data=pd.read_csv(os.path.join(result_path,"substrate_consumption.csv"),header=0,index_col=0)
C=list(no_noise_data["Concentration"])
noise_C={
    "0":C
}
pattern = re.compile(r"\d+.?\d+")
noise_list=[]
for file in os.listdir(result_path):
    noise_str=pattern.findall(file)
    if noise_str:
        noise_list.append(noise_str[0])

print(noise_C)


# 读取噪声数据

for bais in noise_list:
    bias_data=pd.read_csv(os.path.join(result_path,f"substrate_consumption_{bais}.csv"),header=0,index_col=None,usecols=["Concentration"])
    noise_C[f"{bais}"]=list(bias_data["Concentration"])





# 按照积分法给的方法来绘制这个曲线
noise_C_X=defaultdict(list)
noise_C_Y=defaultdict(list)


fig, ax = plt.subplots(5, 2,figsize=(7.7,9.7))
fig.subplots_adjust(0.09, 0.1, 0.99, 0.99)
fig.subplots_adjust(wspace=0.262,hspace=0.364)

lines_list=[]


t=[t*5 for t in range(len(noise_C["0"]))]

for index,(key,value) in enumerate(noise_C.items()):
    for time,C_t in zip(t[1:],value[1:]):
        noise_C_X[f"{key}"].append((value[0] - C_t) / time)
        noise_C_Y[f"{key}"].append(1 / time * log(value[0] / C_t))


ax[0][0].plot(noise_C_X["0"],noise_C_Y["0"],label='A')
ax[0][0].annotate(f"({chr(65)})",xy=(0.9, 0.78),xycoords='axes fraction')

for index,key in enumerate(noise_list):
    X=noise_C_X[f"{key}"]
    Y=noise_C_Y[f"{key}"]
    reg = linear_model.LinearRegression()
    reg.fit(np.array(X).reshape(-1,1),Y)
    print(f"当噪声为{key}时，拟合的斜率和截距分别为:{reg.coef_}和{reg.intercept_}", end=";")
    print(f"计算得k={1/reg.coef_[0]:0.4f},u={1/reg.coef_[0]*reg.intercept_:0.4f}")
    row=int((index+1)/2)
    col=index+1-row*2
    ax[row][col].plot(X,Y,lw=1#,marker='o',ms=2
                      )
    ax[row][col].set_xlabel("X")
    ax[row][col].set_ylabel("Y",labelpad=3)
    ax[row][col].annotate(f"({chr(65+index+1)})", xy=(0.9, 0.78), xycoords='axes fraction')

plt.savefig("the_error_of_the_linear_solution_model.jpg")
plt.show()


