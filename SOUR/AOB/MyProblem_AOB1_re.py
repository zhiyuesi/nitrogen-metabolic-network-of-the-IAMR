# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea
import pandas as pd
"""
相对误差计算AOB的动力学
"""
class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self):
        name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        maxormins = [1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 2 # 初始化Dim（决策变量维数），先是μmax，再是Ks
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [0.4,0.10] # 决策变量下界
        ub = [3.5,1.0] # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        u = Vars[:, [0]]#μmax
        K = Vars[:, [1]]#Ks
        real_data_raw= pd.read_csv('./data/AOB1.csv',usecols=['DO'],encoding='utf-8')

        real_data=real_data_raw[70:126]

        #现从do=4.0 mg/L
        predict_data = []
        for DO in real_data['DO']:
            # 以sec为微元单位了
            for sec in range(5):
                DO = DO - DO * u / (K + DO) / 60
            predict_data.append(DO)

        score=0
        # 因为最后一个预测值是基于最后一个真实数据计算出来的。在real_data中没有对应的值，所以这里要做一个len-1。
        for index in range(len(predict_data) - 1):
            score = score + ((real_data.iat[index + 1, 0] - predict_data[index])/real_data.iat[index + 1, 0]) ** 2 #

        score=score/(len(predict_data) - 1)#REMSE
        pop.ObjV = score # 目标函数值，赋值给pop种群对象的ObjV属性
    