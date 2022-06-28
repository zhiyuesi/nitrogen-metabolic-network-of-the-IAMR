# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea
import pandas as pd

"""
"""
class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        maxormins = [1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 2 # 初始化Dim（决策变量维数），先是μmax，再是Ks
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [0.01,0.00] # 决策变量下界
        ub = [2,10] # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        u = Vars[:, [0]]#μmax
        K = Vars[:, [1]]#Ks
        real_data_raw = pd.read_csv('./data/NOB3.csv',
                                    header=0, index_col=None,
                                    usecols=['Mins', 'NO2--N'],  # 这里要使用两列数据
                                    encoding='utf-8')
        real_data = real_data_raw[1:-1]
        C_s_list = list(real_data["NO2--N"])
        delta_time = []  # 计算每个数据点之前的时间间隔
        for index in real_data.index[:-1]:
            # 时间间隔跟数据点的关系：time(data_i+1)=time(data_i)+deltatime_i
            delta_time.append(real_data.at[index + 1, "Mins"] - real_data.at[index, "Mins"])
        predict_data = []
        # 每个数据点都作为起点进行预测
        for index, C_s in enumerate(C_s_list[:-1]):
            for i in range(delta_time[index]):
                # 这里是每分钟模拟一次
                C_s = C_s - C_s * u / (K + C_s)
            predict_data.append(C_s)
        score = 0
        for index,predict in enumerate(predict_data):# 因为没有使用real_data最后一个数据，所以预测出的predict_data[-1]==C_s_list[-1]
            #score = score + (C_s_list[index + 1] - predict) ** 2  # SSE
            score = score + ((C_s_list[index + 1] - predict)/C_s_list[index + 1]) ** 2
        score=score/(len(predict_data)-1)
        pop.ObjV = score # 目标函数值，赋值给pop种群对象的ObjV属性
    