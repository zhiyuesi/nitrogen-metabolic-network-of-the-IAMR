#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author  : 月司
# @Email   : 815603884@qq.com
# @File    : 上传GitHub\model\simulation_nitrogen_transfer.py
# @version  : 1.0


# 这是以秒为单位来进行模拟的
# 进行物料的计算 dc/dt=input-output-react


import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
from math import fabs
import datetime
import os
import json
import logging
# 并发提高仿真的速度
from concurrent import futures
# 用来获得函数的形参信息
import inspect
import copy

# 日志设置
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
                    )
logger = logging.getLogger(__name__)

# 统计一个周期内的氮底物的变化
ammonia_transfer = {
    "ammonia_oxidation": 0,
    'anammox': 0,
}
nitrite_transfer = {
    "nitrite_oxidation": 0,
    'anammox': 0,
    'nitrite_denitrification': 0,
    "nitrate_respiration": 0,
}
nitrate_transfer = {
    "nitrate_denitrification": 0,
    "nitrate_respiration": 0,
    "nitrite_oxidation": 0,
    'anammox': 0,

}

# 载入动力学数据
kinetics_data = pd.read_json("../kinetics.json", orient="index")

# --------------COD底物在有限时，在硝酸盐反硝化、亚硝酸盐反硝化、硝酸盐呼吸过程的分配关系-------------------#
DE_NITRITE_COD = 1.71
DE_NITRATE_COD = 2.86
NITRATE_RESPIRATION_COD = 1.15
# 载入厌氧的反应速率
filename = "../COD_anaerobic_digestion/COD_anaerobic_digestion_rate.json"
# 载入厌氧COD消化的文件
with open(filename) as fp:
    COD_ANAEROBIC_DIGESTION_RATE = json.load(fp)['COD_anaerobic_digestion_rate']/60
# 硝酸盐呼吸的比例：
filename = "../Anammox&nitrate respiration/mean_nitrate_respiration_proportion.json"
with open(filename, "r", encoding='utf-8') as fp:
    RESPIRATION_PROPORTION = json.load(fp)["mean_nitrate_respiration_proportion"]
COD_OXIDATION_RATE_BY_SOUR = kinetics_data.at["('SOUR', 'CAO')", "umax"]/60
DE_NITRITE_MAX_RATE = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO2')", "umax"]/60
DE_NITRATE_MAX_RATE = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO3')", "umax"] * (
        1 - RESPIRATION_PROPORTION)/60
NITRATE_RESPIRATION_MAX_RATE = kinetics_data.at[
                                   "('Nitrogen_Substance', 'Denitrification_NO3')", "umax"] * RESPIRATION_PROPORTION/60

def COD_FRAC_FOR_DE_NITRITE(DO,VSS,nitrite,nitrate,IsMicroaerobic=True):
    # 判断当前的状态
    if IsMicroaerobic:
        condition="microaerobic"
    else:
        condition="anaerobic"

    de_nitrite_rate=nitrite_denitrification(DO,VSS,nitrite,condition)
    de_nitrate_rate = nitrate_denitrification(DO, VSS, nitrate,RESPIRATION_PROPORTION, condition)
    nitrate_respiration_rate=nitrate_respiration(DO, VSS, nitrate,RESPIRATION_PROPORTION, condition)
    cod_oxidation_rate=COD_oxidation(DO,VSS)
    cod_anaerobic_digestion_rate=COD_ANAEROBIC_DIGESTION_RATE* VSS
    return  de_nitrite_rate * DE_NITRITE_COD / (cod_oxidation_rate*IsMicroaerobic +
                                                                  cod_anaerobic_digestion_rate*(not IsMicroaerobic) +
                                                                  de_nitrate_rate * DE_NITRATE_COD +
                                                                  nitrate_respiration_rate * NITRATE_RESPIRATION_COD +
                                                                  de_nitrite_rate * DE_NITRITE_COD)

def COD_FRAC_FOR_DE_NITRATE(DO,VSS,nitrite,nitrate,IsMicroaerobic=True):
    # 判断当前的状态
    if IsMicroaerobic:
        condition = "microaerobic"
    else:
        condition = "anaerobic"

    de_nitrite_rate = nitrite_denitrification(DO, VSS, nitrite, condition)
    de_nitrate_rate = nitrate_denitrification(DO, VSS, nitrate, RESPIRATION_PROPORTION, condition)
    nitrate_respiration_rate = nitrate_respiration(DO, VSS, nitrate, RESPIRATION_PROPORTION, condition)
    cod_oxidation_rate = COD_oxidation(DO, VSS)
    cod_anaerobic_digestion_rate = COD_ANAEROBIC_DIGESTION_RATE * VSS

    return de_nitrate_rate * DE_NITRATE_COD / (cod_oxidation_rate*IsMicroaerobic +
                                                                  cod_anaerobic_digestion_rate*(not IsMicroaerobic) +
                                                                  de_nitrate_rate * DE_NITRATE_COD +
                                                                  nitrate_respiration_rate * NITRATE_RESPIRATION_COD +
                                                                  de_nitrite_rate * DE_NITRITE_COD)
def COD_FRAC_FOR_DE_NITRATE_RESPIRATION(DO,VSS,nitrite,nitrate,IsMicroaerobic=True):
    # 判断当前的状态
    if IsMicroaerobic:
        condition = "microaerobic"
    else:
        condition = "anaerobic"

    de_nitrite_rate = nitrite_denitrification(DO, VSS, nitrite, condition)
    de_nitrate_rate = nitrate_denitrification(DO, VSS, nitrate, RESPIRATION_PROPORTION, condition)
    nitrate_respiration_rate = nitrate_respiration(DO, VSS, nitrate, RESPIRATION_PROPORTION, condition)
    cod_oxidation_rate = COD_oxidation(DO, VSS)
    cod_anaerobic_digestion_rate = COD_ANAEROBIC_DIGESTION_RATE * VSS
    return nitrate_respiration_rate * NITRATE_RESPIRATION_COD / (cod_oxidation_rate*IsMicroaerobic +
                                                                  cod_anaerobic_digestion_rate*(not IsMicroaerobic) +
                                                                  de_nitrate_rate * DE_NITRATE_COD +
                                                                  nitrate_respiration_rate * NITRATE_RESPIRATION_COD +
                                                                  de_nitrite_rate * DE_NITRITE_COD)

def diluent(substance, diluent_radio=8.0 * 60 * 60):
    '''
    反应器是完全混合的，所以进水进入反应器立刻被稀释
    计算方式：
        Qin(L/h)*HRT(h)=V(L)
        dCin/dt=QCin/V
        联立有：dCin/dt=Cin/HRT(mg·L-1·h-1)

        dCin/dt=substance/diluent_radio
    :param substance:
    :param diluent_radio:HRT(S)
    :return:
    '''
    # 稀释率的计算
    diluent_radio = diluent_radio
    return substance / diluent_radio


def ammonia_oxidation(DO, VSS, ammonia):
    # NH4+-N被AOB氧化的方程,计算消耗的NH4+-N
    mu = kinetics_data.at["('Nitrogen_Substance', 'AOB')", "umax"] * VSS/60

    k = kinetics_data.at["('Nitrogen_Substance', 'AOB')", "ks"]
    k_O = kinetics_data.at["('SOUR', 'AOB')", "ks"]
    # mu_O=kinetics_data.at["('SOUR', 'AOB')","umax"]
    ammonia = ammonia * mu / (k + ammonia) * (DO / (DO + k_O))
    # print("ammonia_oxidation:%.2f"%ammonia)
    return ammonia


def nitrite_oxidation(DO, VSS, nitrite):
    mu = kinetics_data.at["('Nitrogen_Substance', 'NOB')", "umax"] * VSS/60
    k = kinetics_data.at["('Nitrogen_Substance', 'NOB')", "ks"]
    k_O = kinetics_data.at["('SOUR', 'NOB')", "ks"]
    # mu_O=kinetics_data.at["('SOUR', 'NOB')","umax"]
    nitrite = nitrite * mu / (k + nitrite) * (DO / (DO + k_O))
    # print("nitrite_oxidation:%.2f"%nitrite)
    return nitrite


def anammox(VSS, nitrite, inhibition_factor=0.75):
    mu = kinetics_data.at["('Nitrogen_Substance', 'ANAMMOX')", "umax"] * VSS/60
    k = kinetics_data.at["('Nitrogen_Substance', 'ANAMMOX')", "ks"]
    inhibition_factor = inhibition_factor
    nitrite = nitrite * mu / (k + nitrite) * inhibition_factor
    # print("anammox:%.2f"%nitrite)
    return nitrite


def nitrate_denitrification(DO, VSS, nitrate,proportion, conditon='anaerobic'):
    mu = DE_NITRATE_MAX_RATE*VSS

    k = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO3')", "ks"]
    k_O = 0.2
    if conditon == 'anaerobic':
        nitrate = nitrate * mu / (k + nitrate)
    else:
        nitrate = nitrate * mu / (k + nitrate) * (k_O / (k_O + DO))

    # print("nitrate_denitrification:%.2f"%nitrate)
    return nitrate


def nitrate_respiration(DO, VSS, nitrate, proportion, conditon='anaerobic'):
    mu = NITRATE_RESPIRATION_MAX_RATE * VSS
    k = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO3')", "ks"]
    k_O = 0.2
    if conditon == 'anaerobic':
        nitrate = nitrate * mu / (k + nitrate)
    else:
        nitrate = nitrate * mu / (k + nitrate) * (k_O / (k_O + DO))

    # print("nitrate_denitrification:%.2f"%nitrate)
    return nitrate


def nitrite_denitrification(DO, VSS, nitrite, conditon='anaerobic'):
    mu = DE_NITRITE_MAX_RATE * VSS
    k = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO2')", "ks"]
    k_O = 0.2
    if conditon == 'anaerobic':
        nitrite = nitrite * mu / (k + nitrite)
    else:
        nitrite = nitrite * mu / (k + nitrite) * (k_O / (k_O + DO))
    # print("nitrite_denitrification:%.2f"%nitrite)
    return nitrite


def COD_oxidation(DO, VSS, Islimited=False):
    '''
    COD好氧氧化过程
    :param DO:
    :param Islimited: COD底物是否受限
    :return:
    '''
    mu = COD_OXIDATION_RATE_BY_SOUR * VSS
    k = kinetics_data.at["('SOUR', 'CAO')", "ks"]
    if not Islimited:
        COD = DO * mu / (DO + k)

    return COD


def COD_anaerobic_digestion(VSS, COD, Islimited=False):
    mu_constant = COD_ANAEROBIC_DIGESTION_RATE
    if not Islimited:
        COD = mu_constant * VSS
    return COD


def calulate_microaerobic_condition(COD_input, ammonia_input, DO, VSS, COD, ammonia, nitrite, nitrate):
    '''
    :param COD_input:
    :param ammonia_input:
    :param DO:
    :param VSS:
    :param COD:
    :param ammonia:
    :param nitrite:
    :param nitrate:
    :return:
    '''
    COD_current = COD
    COD_aerobic = COD_oxidation(DO, VSS)
    COD_anaerobic = COD_anaerobic_digestion(VSS, COD)
    nitrite_rate = nitrite_denitrification(DO, VSS, nitrite, conditon='microaerobic')
    nitrate_rate = nitrate_denitrification(DO, VSS, nitrate,RESPIRATION_PROPORTION, conditon='microaerobic')
    nitrate_respiration_rate=nitrate_respiration(DO, VSS, nitrate,RESPIRATION_PROPORTION, conditon='microaerobic')


    total_COD_next_reaction = COD_aerobic + COD_anaerobic + nitrate_rate * DE_NITRATE_COD \
                              + nitrate_respiration_rate* NITRATE_RESPIRATION_COD + nitrite_rate * DE_NITRITE_COD

    if COD_current <= total_COD_next_reaction:
        # COD的量不足，是限制性条件
        next_COD = COD_input

        # 计算anammox过程的NO3--N生成
        anammox_nitrate = anammox(VSS, nitrite) * 0.26 / 1.32
        nitrate_transfer['anammox'] += anammox_nitrate
        # 计算亚硝酸盐氧化过程的NO3--N的生成
        NOB_nitrate = nitrite_oxidation(DO, VSS, nitrite)
        nitrate_transfer["nitrite_oxidation"] += NOB_nitrate

        # NO3--N去除
        nitrate2n2 = COD_current * COD_FRAC_FOR_DE_NITRATE(DO,VSS,nitrite,nitrate) / DE_NITRATE_COD
        nitrate2nitrite = COD_current * COD_FRAC_FOR_DE_NITRATE_RESPIRATION(DO,VSS,nitrite,nitrate) / NITRATE_RESPIRATION_COD
        # 硝酸盐反硝化以及硝酸盐呼吸去除硝酸盐的量
        DE_nitrate = nitrate2n2 + nitrate2nitrite
        nitrate_transfer["nitrate_denitrification"] += nitrate2n2
        nitrate_transfer["nitrate_respiration"] += nitrate2nitrite
        # 计算nitrate下一个模拟时刻的浓度
        next_nitrate = nitrate + anammox_nitrate + NOB_nitrate - DE_nitrate - diluent(nitrate)

        # 计算NO2--N生成
        AOB_nitrite = ammonia_oxidation(DO, VSS, ammonia)
        # 计算NO2--N去除
        anammox_nitrite = anammox(VSS, nitrite)
        nitrite_transfer['anammox'] += anammox_nitrite
        NOB_nitrite = nitrite_oxidation(DO, VSS, nitrite)
        nitrite_transfer["nitrite_oxidation"] += NOB_nitrite
        DE_nitrite = COD_current * COD_FRAC_FOR_DE_NITRITE(DO,VSS,nitrite,nitrate) / DE_NITRITE_COD  # 按照比例使用COD
        nitrite_transfer['nitrite_denitrification'] += DE_nitrite
        # 从硝酸盐呼吸来的亚硝酸盐
        nitrite_transfer["nitrate_respiration"] += nitrate2nitrite
        # 计算nitrite下一个模拟时刻的浓度
        next_nitrite = nitrite + AOB_nitrite + nitrate2nitrite - anammox_nitrite - NOB_nitrite - DE_nitrite - diluent(
            nitrite)
    else:
        # COD充足，不是限制性条件
        next_COD = COD_input + COD - total_COD_next_reaction

        # 更新硝酸盐
        anammox_nitrate = anammox(VSS, nitrite) * 0.26 / 1.32
        nitrate_transfer['anammox'] += anammox_nitrate
        NOB_nitrate = nitrite_oxidation(DO, VSS, nitrite)
        nitrate_transfer["nitrite_oxidation"] += NOB_nitrate

        nitrate2n2 = nitrate_rate
        nitrate2nitrite = nitrate_respiration_rate
        DE_nitrate = nitrate2nitrite + nitrate2n2
        nitrate_transfer["nitrate_denitrification"] += nitrate2n2
        nitrate_transfer["nitrate_respiration"] += nitrate2nitrite
        next_nitrate = nitrate + anammox_nitrate + NOB_nitrate - DE_nitrate - diluent(nitrate)

        # 更新亚硝酸盐
        AOB_nitrite = ammonia_oxidation(DO, VSS, ammonia)
        anammox_nitrite = anammox(VSS, nitrite)
        nitrite_transfer['anammox'] += anammox_nitrite
        NOB_nitrite = nitrite_oxidation(DO, VSS, nitrite)
        nitrite_transfer["nitrite_oxidation"] += NOB_nitrite
        DE_nitrite = nitrite_rate
        nitrite_transfer['nitrite_denitrification'] += DE_nitrite
        nitrite_transfer["nitrate_respiration"] += nitrate2nitrite
        next_nitrite = nitrite + AOB_nitrite + nitrate2nitrite - anammox_nitrite - NOB_nitrite - DE_nitrite - diluent(
            nitrite)

    # 最后更新NH4+-N的浓度变化
    anammox_ammonia = anammox(VSS, nitrite) / 1.32
    AOB_ammonia = ammonia_oxidation(DO, VSS, ammonia)
    ammonia_transfer['anammox'] += anammox_ammonia
    ammonia_transfer["ammonia_oxidation"] += AOB_ammonia
    next_ammonia = ammonia + ammonia_input - diluent(ammonia) - AOB_ammonia - anammox_ammonia
    return next_COD, next_nitrite, next_nitrate, next_ammonia


def calculate_anaerobic_condition(COD_input, ammonia_input, DO, VSS, COD, ammonia, nitrite, nitrate):
    '''

    :param COD_input:
    :param ammonia_input:
    :param DO:
    :param VSS:
    :param COD:
    :param ammonia:
    :param nitrite:
    :param nitrate:
    :return:
    '''
    COD_current = COD

    COD_anaerobic = COD_anaerobic_digestion(VSS, COD)
    nitrite_rate = nitrite_denitrification(DO, VSS, nitrite)
    nitrate_rate = nitrate_denitrification(DO, VSS, nitrate,RESPIRATION_PROPORTION)
    nitrate_respiration_rate=nitrate_respiration(DO, VSS, nitrate,RESPIRATION_PROPORTION)

    total_COD_next_reaction = COD_anaerobic + nitrate_rate * DE_NITRATE_COD \
                                  + nitrate_respiration_rate* NITRATE_RESPIRATION_COD + nitrite_rate * DE_NITRITE_COD

    if COD_current <= total_COD_next_reaction:
        # COD的量不足，是限制性条件
        next_COD = COD_input

        anammox_nitrate = anammox(VSS, nitrite, inhibition_factor=1.0) * 0.26 / 1.32
        nitrate_transfer['anammox'] += anammox_nitrate


        nitrate2n2 = COD_current * COD_FRAC_FOR_DE_NITRATE(DO,VSS,nitrite,nitrate,IsMicroaerobic=False) / DE_NITRATE_COD
        nitrate2nitrite = COD_current * COD_FRAC_FOR_DE_NITRATE_RESPIRATION(DO,VSS,nitrite,nitrate,IsMicroaerobic=False) / NITRATE_RESPIRATION_COD

        DE_nitrate = nitrate2n2 + nitrate2nitrite
        nitrate_transfer["nitrate_denitrification"] += nitrate2n2
        nitrate_transfer["nitrate_respiration"] += nitrate2nitrite

        next_nitrate = nitrate + anammox_nitrate - DE_nitrate - diluent(nitrate)


        anammox_nitrite = anammox(VSS, nitrite, inhibition_factor=1.0)
        nitrite_transfer['anammox'] += anammox_nitrite
        DE_nitrite = COD_current * COD_FRAC_FOR_DE_NITRITE(DO,VSS,nitrite,nitrate,IsMicroaerobic=False) / DE_NITRITE_COD
        nitrite_transfer['nitrite_denitrification'] += DE_nitrite
        nitrite_transfer["nitrate_respiration"] += nitrate2nitrite
        next_nitrite = nitrite + nitrate2nitrite - anammox_nitrite - DE_nitrite - diluent(nitrite)
    else:
        # COD充足，不是限制性条件
        next_COD = COD_input + COD - total_COD_next_reaction
        # 更新硝酸盐
        anammox_nitrate = anammox(VSS, nitrite, inhibition_factor=1.0) * 0.26 / 1.32
        nitrate_transfer['anammox'] += anammox_nitrate

        nitrate2n2 = nitrate_rate * (1 - RESPIRATION_PROPORTION)
        nitrate2nitrite = nitrate_rate * RESPIRATION_PROPORTION
        DE_nitrate = nitrate2nitrite + nitrate2n2
        nitrate_transfer["nitrate_denitrification"] += nitrate2n2
        nitrate_transfer["nitrate_respiration"] += nitrate2nitrite
        next_nitrate = nitrate + anammox_nitrate - DE_nitrate - diluent(nitrate)

        # 更新亚硝酸盐
        anammox_nitrite = anammox(VSS, nitrite, inhibition_factor=1.0)
        nitrite_transfer['anammox'] += anammox_nitrite
        DE_nitrite = nitrite_rate
        nitrite_transfer['nitrite_denitrification'] += DE_nitrite
        nitrite_transfer["nitrate_respiration"] += nitrate2nitrite
        next_nitrite = nitrite + nitrate2nitrite - anammox_nitrite - DE_nitrite - diluent(nitrite)

    anammox_ammonia = anammox_nitrite / 1.32
    ammonia_transfer['anammox'] += anammox_ammonia
    next_ammonia = ammonia + ammonia_input - diluent(ammonia) - anammox_ammonia
    return next_COD, next_nitrite, next_nitrate, next_ammonia


def simulate(ammonia_input, COD_input,
             ammonia_init, nitrite_init, nitrate_init, COD_init,
             VSS=5.5, DO=0.2,
             anaerobic_period=6, microaerobic_period=18,
             result_path="", TrackList=None,
             accuracy=0.01, figure=False) -> list:
    '''
    :param ammonia_input:
    :param COD_input:
    :param ammonia_init:
    :param nitrite_init:
    :param nitrate_init:
    :param COD_init:
    :param VSS:
    :param DO:
    :param anaerobic_period:
    :param microaerobic_period:
    :param result_path:
    :param TrackList:要追踪底物变化的列表
    :param accuracy:
    :param figure:
    :return result_list: 返回的是输入的NH4+-N以及COD浓度，以及稳定周期时TrackList中设定物质的出水浓度(TrackList默认为[NH4+-N、NO2--N、NO3--N]）
    '''

    # 跟踪底物的浓度随时间的变化，默认只跟踪三氮的变化
    if not TrackList:
        TrackList = ['ammonia', 'nitrite', 'nitrate']
    # 设为带默认值的字典类型数据
    CurrentPeriod_data = defaultdict(list)
    LastPeriod_data = defaultdict(list)
    CurrentPeriod_Mean = defaultdict(float)
    LastPeriod_Mean = defaultdict(float)
    #
    if not result_path:
        result_path = './result'
    # 当不存在result路径，需要提前建立
    if not os.path.exists(result_path):
        # 由于并发，极小概论会两个进程在极短的时间内，同时进行目录创建，会发生错误
        try:
            os.makedirs(result_path)
        except FileExistsError:
            pass
    # 过程仿真数据输出到csv文件中，文件的命名规则
    calulate_trace = os.path.join(result_path,
                                  f'{COD_input:.2f}_{ammonia_input:.2f}_simulation_data_{datetime.datetime.now().strftime("%Y-%m-%d %H%M%S")}.csv')
    with open(calulate_trace, 'w', encoding='utf-8') as f:
        f.write(','.join(['time', 'ammonia', 'nitrite', 'nitrate', 'COD']))
        f.write('\n')
    # 将每个仿真，最后周期的各个功能菌的脱氮贡献输出到文件中。
    removal_trace = os.path.join(result_path,
                                 f'{COD_input:.2f}_{ammonia_input:.2f}_removal_data_{datetime.datetime.now().strftime("%Y-%m-%d %H%M%S")}.json')

    ammonia_input = diluent(ammonia_input)  # 计算出每秒的进水率Ds/Dt
    COD_input = diluent(COD_input)
    print(ammonia_input, COD_input)

    # 这里是一个曝气周期的情况。1代表曝气状态，0代表不曝气状态
    simulation_period = [1] * microaerobic_period + [0] * anaerobic_period
    # 就是显示时间列表,用来画图
    show_time = [0]
    # 代表了已经运行的分钟数
    time_has_gone = 0
    # 下面是反应器中初始值的浓度
    COD = COD_init
    ammonia = ammonia_init
    nitrite = nitrite_init
    nitrate = nitrate_init

    if figure:
        fig, ax = plt.subplots(4, 1, sharex=True)
        fig.subplots_adjust(hspace=0)
        plt.ion()
        ammonia_list = [ammonia]
        nitrite_list = [nitrite]
        nitrate_list = [nitrate]
        COD_list = [COD]

        line_ammonia, = ax[0].plot(show_time, ammonia_list)
        line_nitrite, = ax[1].plot(show_time, nitrite_list)
        line_nitrate, = ax[2].plot(show_time, nitrate_list)
        line_COD, = ax[3].plot(show_time, COD_list)
    while True:
        # 按周期计算污染物的降解规律，最终是按秒计算的
        for time in simulation_period:

            # 以1min为步长更新数据
            if time == 1:  # 微氧状态
                for sec in range(0, 60):
                    COD, nitrite, nitrate, ammonia = calulate_microaerobic_condition(COD_input, ammonia_input, DO, VSS,
                                                                                     COD,
                                                                                     ammonia, nitrite, nitrate)
            else:
                for sec in range(0, 60):  # 厌氧状态

                    COD, nitrite, nitrate, ammonia = calculate_anaerobic_condition(COD_input, ammonia_input, 0, VSS,
                                                                                   COD,
                                                                                   ammonia, nitrite, nitrate)

            print("{0:.2f}  {1:.2f}  {2:.2f}  {3:.2f}".format(ammonia, nitrite, nitrate, COD))

            time_has_gone += 1
            if figure:
                show_time.append(time_has_gone)
                ammonia_list.append(ammonia)
                nitrite_list.append(nitrite)
                nitrate_list.append(nitrate)
                COD_list.append(COD)
                # 更新数据图
                line_ammonia.set_xdata(show_time)
                line_ammonia.set_ydata(ammonia_list)
                line_nitrite.set_xdata(show_time)
                line_nitrite.set_ydata(nitrite_list)
                line_nitrate.set_xdata(show_time)
                line_nitrate.set_ydata(nitrate_list)
                line_COD.set_xdata(show_time)
                line_COD.set_ydata(COD_list)
                ax[0].set_ylim(-0.5, 50)
                ax[1].set_ylim(-0.05, 2)
                ax[2].set_ylim(0, 50)
                ax[3].set_ylim(0, 10)
                ax[0].set_xlim(0, len(show_time) + 10)
                plt.pause(0.001)
            # 保存一个周期内的数据
            for substrate in TrackList:
                CurrentPeriod_data[substrate].append(eval(substrate))

            with open(calulate_trace, 'a', encoding='utf-8') as f:
                f.write(','.join([str(x) for x in [time_has_gone, ammonia, nitrite, nitrate, COD]]) + '\n')

        print(
            "ammonia removal:\nammonia oxidation:{0:.4f}\tanammox:{1:.4f}".format(ammonia_transfer["ammonia_oxidation"],
                                                                                  ammonia_transfer['anammox']))
        print("nitrite removal:\nnitrite_oxidation:{0:.4f}\tanammox:{1:.4f}\tnitrite_denitrification{2:.4f}".format(
            nitrite_transfer["nitrite_oxidation"], nitrite_transfer['anammox'],
            nitrite_transfer['nitrite_denitrification']))
        print("nitrate removal:\nnitrate_denitrification{0:.4f}".format(nitrate_transfer['nitrate_denitrification']))
        print("nitrate generate:\nnitrate_oxidation{0:.2f}\tanammox:{1:.02f}".format(
            nitrate_transfer["nitrite_oxidation"], nitrate_transfer['anammox']))

        for substrate in TrackList:
            CurrentPeriod_Mean[substrate] = np.mean(CurrentPeriod_data[substrate])
        # 判断是否反应稳定
        if not LastPeriod_Mean:
            LastPeriod_Mean = CurrentPeriod_Mean.copy()
        else:
            # 判断是否都满足了条件
            result_flags = [
                fabs(CurrentPeriod_Mean[substrate] - LastPeriod_Mean[substrate]) <= accuracy for
                substrate in TrackList]
            # 当仿真时间小于14小时，且存在污染物一个指标误差大于0.001，则再次计算
            if False in result_flags and time_has_gone < 60 * 14:
                for key in ammonia_transfer:
                    ammonia_transfer[key] = 0
                for key in nitrite_transfer:
                    nitrite_transfer[key] = 0
                for key in nitrate_transfer:
                    nitrate_transfer[key] = 0

                LastPeriod_Mean = CurrentPeriod_Mean.copy()
                CurrentPeriod_data.clear()

            else:
                with open(removal_trace, 'w', encoding="utf-8") as file_obj:
                    json.dump([ammonia_transfer, nitrite_transfer, nitrate_transfer], file_obj,
                              ensure_ascii=False)
                plt.close()
                result_list=[ammonia_input,COD_input]+[CurrentPeriod_Mean[substrate] for substrate in TrackList]
                return result_list

def simulate_proxy(params:dict):
    # 获得其形参
    target_args=inspect.signature(simulate).parameters
    target_args_set=set(target_args)
    params_keys_set=set(params.keys())
    # 将params里面，目标函数不存在的函数给剔除
    unknown_params=params_keys_set-target_args_set
    for key in unknown_params:
        del params[key]

    # 将目标函数存在的形参，但是输入的参数并没有包括的参数，即使用关键字参数的默认值
    default_params=target_args_set-params_keys_set
    for key in default_params:
        params[key]=target_args[key].default
    #print(params)
    return simulate(ammonia_init=params["ammonia_init"],
                   nitrite_init=params["nitrite_init"],
                   nitrate_init=params["nitrate_init"],
                   COD_init=params["COD_init"],
                   COD_input=params["COD_input"],
                   ammonia_input=params["ammonia_input"],
                   DO=params["DO"],
                    VSS=params["VSS"],
                    figure=params["figure"]
                    )


def simulate_params(experimental_datas:pd.DataFrame,**kwargs):
    '''

    :param experimental_datas:
    :param kwargs:simulate（）其他的关键字参数，必须是list类型，且长度等于experimental_datas.shape[0]
    :return:
    '''
    data_length=experimental_datas.shape[0]
    # print(kwargs)
    for key,value in kwargs.items():
        assert isinstance(value,list),f"{key}不为列表类型数据"
        assert len(value)==data_length,f"{key}长度跟数据记录的长度不符"

    result_list=[]
    # 因为首个数据仅用于确定反应器的初始浓度，所以这里从1这个索引开始循环
    for pos,daily_index in enumerate(experimental_datas.index[1:]):
        temp_dict={
            "ammonia_init": experimental_datas.at[daily_index - 1, 'Ammonia_OUT'],
            "nitrite_init": experimental_datas.at[daily_index - 1, 'Nitrite_OUT'],
            "nitrate_init": experimental_datas.at[daily_index - 1, 'Nitrate_OUT'],
            "COD_init": experimental_datas.at[daily_index - 1, 'COD_OUT'],
            "COD_input": experimental_datas.at[daily_index, 'COD_removal'],
            "ammonia_input": experimental_datas.at[daily_index, 'TKN'],
        }
        # 因为enumerate，是从0开始的，但是我们DO和VSS不需要初始化，数据用的是[1:-1]
        for key, value in kwargs.items():
            temp_dict[key]=value[pos+1]

        result_list.append(temp_dict)

    return result_list

def linear_VSS(VSS_datas:pd.Series):
    '''
    通过线性插值的方法计算每日的VSS数据
    :param VSS_datas:
    :return:
    '''
    index_list=list(VSS_datas.index)
    VSS=[]
    for pos in range(len(index_list)-1):
        bin=(VSS_datas[index_list[pos+1]]-VSS_datas[index_list[pos]])/(index_list[pos+1]-index_list[pos])
        for index in range(index_list[pos],index_list[pos+1]):
            ins_index=index-index_list[pos]
            VSS.append(round(ins_index*bin + VSS_datas[index_list[pos]],2))
    # 需要加入最后一个VSS的点
    VSS.append(VSS_datas.iat[-1])
    return VSS


if __name__ == '__main__':

    experimental_data = pd.read_csv("./daily_pollutants.csv", header=0, index_col=False)

    experimental_data['TKN'] = None
    for index in experimental_data.index:

        if experimental_data.at[index, "TN_IN"] > experimental_data.at[index, 'TIN_IN']:
            TON_in = experimental_data.at[index, "TN_IN"] - experimental_data.at[index, 'TIN_IN']
        else:
            TON_in = 0

        if experimental_data.at[index, "TN_OUT"] > experimental_data.at[index, 'TIN_OUT']:
            TON_out = experimental_data.at[index, "TN_OUT"] - experimental_data.at[index, 'TIN_OUT']
        else:
            TON_out = 0
        experimental_data.at[index, 'TKN'] = experimental_data.at[index, 'Ammonia_IN'] + TON_in - TON_out


    experimental_data['COD_removal'] = experimental_data['COD_IN'] - experimental_data['COD_F']

    # 并发线程
    simulate_datas=experimental_data.loc[42:,['Ammonia_OUT','Nitrite_OUT','Nitrate_OUT','COD_OUT','COD_removal','TKN']]
    DO=list(experimental_data.loc[42:,"DO"])
    VSS_datas=linear_VSS(experimental_data["VSS"].dropna().iloc[-3:])
    simulate_params_list=simulate_params(simulate_datas,DO=DO,VSS=VSS_datas)
    # 4核运行
    with futures.ProcessPoolExecutor(max_workers=4) as executor:
        # map是按照迭代的顺序，返回结果。
        for result in executor.map(simulate_proxy,simulate_params_list):
            # 一次保存NH4+-N_input,cod_input,NH4+-N、NO2--N、NO3--N的平均值信息
            with open('./result/simulation_data.csv', mode='a', encoding='utf-8') as f:
                f.write(','.join([str(x) for x in [*result]]) + '\n')

    # 这里是要绘制最后的仿真图
    fig, ax = plt.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    simulation_data = pd.read_csv('./result/simulation_data.csv', encoding='utf-8', header=None)
    simulation_data.columns = ["COD_Input", "Ammonia_Input", "Ammonia_OUT", "Nitrite_OUT", "Nitrate_OUT"]
    show_list = ["Ammonia_OUT", "Nitrite_OUT", "Nitrate_OUT"]
    for (index, substrate) in enumerate(show_list):
        # F代表真实的数据
        experimental_data_label = substrate.split("_")[0] + "_F"
        ax[index].plot(list(experimental_data[experimental_data_label])[-14:], marker='o', color='red',
                       label=experimental_data_label)
        ax[index].plot(list(simulation_data[substrate])[-14:], color='blue', marker='o', label=substrate)
        ax[index].legend()
    plt.legend()
    plt.show()
