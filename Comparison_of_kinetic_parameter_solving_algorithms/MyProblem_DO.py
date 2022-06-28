# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea
import pandas as pd
"""

"""
class MyProblem(ea.Problem): # 继承Problem父类
    def __init__(self,datafile):
        self.datafile=datafile
        print("\t"*4+self.datafile)
        name = 'MyProblem' # 初始化name（函数名称，可以随意设置）
        M = 1 # 初始化M（目标维数）
        maxormins = [1] # 初始化maxormins（目标最小最大化标记列表，1：最小化该目标；-1：最大化该目标）
        Dim = 2 # 初始化Dim（决策变量维数），先是μmax，umax的单位是mg/L/s;再是Ks
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [0.01,0.10] # 决策变量下界
        ub = [0.4,1.0] # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        u = Vars[:, [0]]#μmax
        K = Vars[:, [1]]#Ks
        real_data = pd.read_csv(self.datafile, delimiter=',', encoding='utf-8', usecols=['Concentration'])

        #从do=4.0 mg/L开始
        predict_data = []
        for DO in real_data['Concentration']:
            #以sec为微元单位了
            for sec in range(5):
                DO = DO - DO * u / (K + DO)
            predict_data.append(DO)
        score=0

        for index in range(len(predict_data) - 1):
            # 最后一个预测数据，是用真实数据最后一个数据点预测得到的，所以在真实数据列表中没有对应值，所以这里取的是len-1
            score = score + ((real_data.iat[index + 1, 0] - predict_data[index])/real_data.iat[index + 1, 0]) ** 2  # SSE
        score=score/(len(predict_data) - 2)#使其变化比较显著
        pop.ObjV = score # 目标函数值，赋值给pop种群对象的ObjV属性
    