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
        Dim = 1 #这里只有一个自变量，那就是硝酸盐呼吸的比例。
        varTypes = [0] * Dim # 初始化varTypes（决策变量的类型，元素为0表示对应的变量是连续的；1表示是离散的）
        lb = [0.0,] # 决策变量下界
        ub = [1.0,] # 决策变量上界
        lbin = [1] * Dim # 决策变量下边界（0表示不包含该变量的下边界，1表示包含）
        ubin = [1] * Dim # 决策变量上边界（0表示不包含该变量的上边界，1表示包含）
        # 调用父类构造方法完成实例化
        ea.Problem.__init__(self, name, M, maxormins, Dim, varTypes, lb, ub, lbin, ubin)
    
    def aimFunc(self, pop): # 目标函数
        Vars = pop.Phen # 得到决策变量矩阵
        f = Vars[:, [0]]#
        # 载入VSS的数据
        VSS_data = pd.read_csv('../VSS.csv', header=0, index_col=[0, 1, 2])
        # 生物量等参数
        VSS = VSS_data.at[("Composite_Processes","Anammox&nitrate_respiration","A&nitrate_respiration1"),"VSS"]
        #载入动力学数据
        kinetics_data = pd.read_json("../kinetics.json", orient="index")

        d_n3_mu = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO3')","umax"]  * VSS
        d_n3_k = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO3')","ks"]

        d_n2_mu = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO2')","umax"]* VSS
        d_n2_k = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO2')","ks"]

        a_n2_mu = kinetics_data.at["('Nitrogen_Substance', 'ANAMMOX')","umax"] * VSS
        a_n2_k = kinetics_data.at["('Nitrogen_Substance', 'ANAMMOX')","ks"]

        real_data_raw = pd.read_csv("./data/A&nitrate_respiration1.csv",header=0,index_col=None,
                                      usecols=["Mins","NH4+-N","NO2--N","NO3--N","COD"])

        real_data = real_data_raw[1:-2]

        # 真实反应物的值，计算score
        C_s_list = list(real_data["NO3--N"])

        real_nh_list=list(real_data["NH4+-N"])
        delta_time = []

        for index in real_data.index[1:]:
            delta_time.append(real_data.at[index, "Mins"] - real_data.at[index - 1, "Mins"])
        #print(delta_time)
        cal_data = []

        # 以第一个点进行预测
        C_s = C_s_list[0]
        # 模拟前的亚硝盐浓度
        no2 = real_data.at[1,"NO2--N"]
        no2_list = [no2, ]
        nh4 = real_data.at[1, "NH4+-N"]
        nh4_list = [nh4, ]
        for delta in delta_time:
            for i in range(delta):
                # 这里是每分钟模拟一次

                r_no3 = C_s * d_n3_mu / (d_n3_k + C_s)
                g_no2 = C_s * d_n3_mu / (d_n3_k + C_s) * f
                r_no2_by_anammox = a_n2_mu * no2 / (a_n2_k + no2)
                r_nh_by_anammox = r_no2_by_anammox / 1.32
                r_no2_by_de = d_n2_mu * no2 / (d_n2_k + no2)
                nh4 = nh4 - r_nh_by_anammox
                g_no3 = r_no2_by_anammox / 1.32 * 0.26
                no2 = no2 + g_no2 - r_no2_by_anammox-r_no2_by_de
                C_s = C_s - r_no3 + g_no3
            nh4_list.append(nh4)
            no2_list.append(no2)
            cal_data.append(C_s)

        score = 0
        for index in range(len(nh4_list)):  # 因为使用real_data最后一个数据，预测出的cal_data，超出了real_data的范围，所以这里要减一
            score = score + ((real_nh_list[index] - nh4_list[index])/real_nh_list[index]) ** 2

        score=score/(len(nh4_list)-1)
        pop.ObjV = score # 目标函数值，赋值给pop种群对象的ObjV属性
    