# -*- coding: utf-8 -*-
'''
@Author  :  月司
@Email   :  815603884@qq.com
@File    :  上传GitHub\Comparison_of_dynamic_parameter_solving_algorithms\Generate_test_data_with_Gaussian_error.py
@version :  1.0
@Desc    : 给标准的莫诺特方程来添加高斯噪声，然后再进行模拟，保存这些数据，绘制了在不同误差下的动图。
'''
import numpy as np
from math import log
import pandas as pd
import os
import matplotlib.pyplot as plt

#用二分法求出每5秒的浓度变化
def target(C, t, C0):
    '''

    :param C:
    :param t:
    :param C0:
    :return:
    '''
    V_max = 0.02  # mg/L/s
    K = 0.36
    return -t + K / V_max * log(C0 / C) + (C0 - C) / V_max #=0


# 使用二分法求极值
def bisection(left_end, right_end, t, C0, epsilon=1e-8):
    '''
    这里使用二分法求方程的零点。
    :param left_end:
    :param right_end:
    :param t: t=0时，C为C0。这里一直以初始时刻为计算点
    :param C0: 初始值，这里一直以初始浓度为计算点
    :param epsilon: 精度
    :return:
    '''

    sign_target_a = np.sign(target(left_end, t, C0))# np.sign只专注符号
    sign_target_b = np.sign(target(right_end, t, C0))
    if sign_target_a==sign_target_b:
        assert "该区间不存在根或偶数根"
    while np.abs(right_end - left_end) > epsilon:
        mid = (left_end + right_end) / 2
        if np.sign(target(mid, t, C0)) == sign_target_a:
            left_end = mid
        else:
            right_end = mid
    return left_end


t = np.arange(0, 305, 5)
C = [4]
for time in t[1:]:
    C0 = 4.0
    print(time)
    # f 计算方程组的误差，[ 1=]是未知数的初始值
    result = bisection(1e-16, C0, time, C0)
    C.append(result)

# 先把模拟数据给保存了

result_data=pd.DataFrame([[time,c_t] for time,c_t in zip(t,C)],columns=["Time","Concentration"])

try:
    result_path = "./data"
    result_data.to_csv(os.path.join(result_path,"substrate_consumption.csv"))
except FileNotFoundError:
    os.mkdir(result_path)
    result_data.to_csv(os.path.join(result_path, "substrate_consumption.csv"))

noise_C={
    "0.0001":[4],
    "0.0003":[4],
    "0.0005":[4],
    "0.001":[4],
    "0.003":[4],
    "0.005":[4],
    "0.01":[4],
    "0.03":[4],
    "0.05":[4],
}

# 给C添加高斯噪声。按照变异系数为5%，计算每个浓度点对应的标准差，然后随机生成噪声。
 # mean and standard deviation

for bais in list(noise_C.keys()):
    for C_t in C[1:]:
        mu, sigma = C_t, C_t * float(bais)
        s = np.random.normal(mu, sigma)
        noise_C[f"{bais}"].append(s)
        result_data = pd.DataFrame([[time, c_t] for time, c_t in zip(t, noise_C[f"{bais}"])], columns=["Time", "Concentration"])
        result_data.to_csv(os.path.join(result_path,f"substrate_consumption_{bais}.csv"))

