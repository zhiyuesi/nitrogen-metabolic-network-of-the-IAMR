# -*- coding: utf-8 -*-
'''
@Author  :  月司
@Email   :  815603884@qq.com
@File    :  matplotlib_config.py
@version :  1.0
@Desc    : matplot绘图设定
'''
import matplotlib.pyplot as plt

plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['font.family'] = ['serif']
plt.rcParams['font.serif']=['Times New Roman']
plt.rcParams['axes.unicode_minus']=False
Legend_Dict={
    "fancybox":False,
    "shadow":False,
    "frameon":False
}