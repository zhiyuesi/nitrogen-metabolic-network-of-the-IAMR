# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea
import pandas as pd
"""
相对误差计算NOB的动力学参数
"""
class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        maxormins = [1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 2 # 初始化Dim（决策变量维数），先是μmax，再是Ks
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [0.01,-0.1] # 决策变量下界
        ub = [0.4,1.2] # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        u = Vars[:, [0]]#μmax
        K = Vars[:, [1]]#Ks
        real_data_raw = pd.read_csv('./data/NOB1.csv', delimiter=',', encoding='utf-8', usecols=['DO'])
        real_data_raw = real_data_raw[258:665]
        real_data = []
        for num in range(0, int(len(real_data_raw['DO']) / 6) + 1):
            real_data.append(real_data_raw.iat[num * 6, 0])
        # 现从do=4.0 mg/L开始
        predict_data = []
        for DO in real_data:
            for i in range(6):
                DO = DO - DO * u / (K + DO) * 5 / 60
            predict_data.append(DO)
        score = 0
        for index, perdict in enumerate(predict_data[:-1]):
            # real_data[index + 1]是预测值所对应的真实值
            score = score + ((real_data[index + 1] - perdict)/real_data[index + 1]) ** 2

        score = score / (len(predict_data) - 1)  # REMSE
        pop.ObjV = score # 目标函数值，赋值给pop种群对象的ObjV属性
    