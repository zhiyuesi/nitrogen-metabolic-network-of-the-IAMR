#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author  : 月司
# @Email   : 815603884@qq.com
# @File    : nitrogen_transfer_no_anaerobic_period.py
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
#日志模块
import logging
# 进度条模块
from tqdm import tqdm
# 命令行模式
import argparse
# 并发提高仿真的速度
from concurrent import futures
# 用来获得函数的形参信息
import inspect
# 正则表达式
import re
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



filename = "../Anammox&nitrate respiration/mean_nitrate_respiration_proportion.json"
with open(filename, "r", encoding='utf-8') as fp:
    RESPIRATION_PROPORTION = json.load(fp)["mean_nitrate_respiration_proportion"]
COD_OXIDATION_RATE_BY_SOUR = kinetics_data.at["('SOUR', 'CAO')", "umax"]/60
DE_NITRITE_MAX_RATE = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO2')", "umax"]/60
DE_NITRATE_MAX_RATE = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO3')", "umax"] * (
        1 - RESPIRATION_PROPORTION)/60
NITRATE_RESPIRATION_MAX_RATE = kinetics_data.at[
                                   "('Nitrogen_Substance', 'Denitrification_NO3')", "umax"] * RESPIRATION_PROPORTION/60

AMMONIA_OXIDATION_MU=kinetics_data.at["('Nitrogen_Substance', 'AOB')", "umax"]
AMMONIA_OXIDATION_K=kinetics_data.at["('Nitrogen_Substance', 'AOB')", "ks"]
AMMONIA_OXIDATION_K_O = kinetics_data.at["('SOUR', 'AOB')", "ks"]
NITRITE_OXIDATION_MU = kinetics_data.at["('Nitrogen_Substance', 'NOB')", "umax"]
NITRITE_OXIDATION_K = kinetics_data.at["('Nitrogen_Substance', 'NOB')", "ks"]
NITRITE_OXIDATION_K_O = kinetics_data.at["('SOUR', 'NOB')", "ks"]
ANAMMOX_MU = kinetics_data.at["('Nitrogen_Substance', 'ANAMMOX')", "umax"]
ANAMMOX_K = kinetics_data.at["('Nitrogen_Substance', 'ANAMMOX')", "ks"]
NITRATE_DENITRIFICATION_K = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO3')", "ks"]
NITRATE_RESPIRATION_K = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO3')", "ks"]
NITRITE_DENITRIFICATION_K = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO2')", "ks"]
COD_OXIDATION_K = kinetics_data.at["('SOUR', 'CAO')", "ks"]
def COD_FRAC_FOR_DE_NITRITE(DO,VSS,nitrite,nitrate,IsMicroaerobic=True,
                            COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR=1.0):
    # 判断当前的状态
    if IsMicroaerobic:
        condition="microaerobic"
    else:
        condition="anaerobic"

    de_nitrite_rate=nitrite_denitrification(DO,VSS,nitrite,condition)
    de_nitrate_rate = nitrate_denitrification(DO, VSS, nitrate,RESPIRATION_PROPORTION, condition)
    nitrate_respiration_rate=nitrate_respiration(DO, VSS, nitrate,RESPIRATION_PROPORTION, condition)
    cod_oxidation_rate=COD_oxidation(DO,VSS)

    cod_anaerobic_digestion_rate = COD_ANAEROBIC_DIGESTION_RATE * VSS*COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR
    return  de_nitrite_rate * DE_NITRITE_COD / (cod_oxidation_rate*IsMicroaerobic +
                                                                  cod_anaerobic_digestion_rate +
                                                                  de_nitrate_rate * DE_NITRATE_COD +
                                                                  nitrate_respiration_rate * NITRATE_RESPIRATION_COD +
                                                                  de_nitrite_rate * DE_NITRITE_COD)

def COD_FRAC_FOR_DE_NITRATE(DO,VSS,nitrite,nitrate,IsMicroaerobic=True,
                            COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR=1.0):
    # 判断当前的状态
    if IsMicroaerobic:
        condition = "microaerobic"
    else:
        condition = "anaerobic"


    de_nitrite_rate = nitrite_denitrification(DO, VSS, nitrite, condition)
    de_nitrate_rate = nitrate_denitrification(DO, VSS, nitrate, RESPIRATION_PROPORTION, condition)
    nitrate_respiration_rate = nitrate_respiration(DO, VSS, nitrate, RESPIRATION_PROPORTION, condition)
    cod_oxidation_rate = COD_oxidation(DO, VSS)

    cod_anaerobic_digestion_rate = COD_ANAEROBIC_DIGESTION_RATE * VSS*COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR

    return de_nitrate_rate * DE_NITRATE_COD / (cod_oxidation_rate*IsMicroaerobic +
                                                                cod_anaerobic_digestion_rate +
                                                                  de_nitrate_rate * DE_NITRATE_COD +
                                                                  nitrate_respiration_rate * NITRATE_RESPIRATION_COD +
                                                                  de_nitrite_rate * DE_NITRITE_COD)
def COD_FRAC_FOR_DE_NITRATE_RESPIRATION(DO,VSS,nitrite,nitrate,IsMicroaerobic=True,
                            COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR=1.0):
    # 判断当前的状态
    if IsMicroaerobic:
        condition = "microaerobic"
    else:
        condition = "anaerobic"

    de_nitrite_rate = nitrite_denitrification(DO, VSS, nitrite, condition)
    de_nitrate_rate = nitrate_denitrification(DO, VSS, nitrate, RESPIRATION_PROPORTION, condition)
    nitrate_respiration_rate = nitrate_respiration(DO, VSS, nitrate, RESPIRATION_PROPORTION, condition)
    cod_oxidation_rate = COD_oxidation(DO, VSS)

    cod_anaerobic_digestion_rate = COD_ANAEROBIC_DIGESTION_RATE * VSS*COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR
    return nitrate_respiration_rate * NITRATE_RESPIRATION_COD / (cod_oxidation_rate*IsMicroaerobic +cod_anaerobic_digestion_rate +
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
    :param diluent_radio:HRT（s）
    :return:
    '''

    diluent_radio = diluent_radio
    return substance / diluent_radio


def ammonia_oxidation(DO, VSS, ammonia):

    mu = AMMONIA_OXIDATION_MU * VSS/60
    k = AMMONIA_OXIDATION_K
    k_O = AMMONIA_OXIDATION_K_O
    # mu_O=kinetics_data.at["('SOUR', 'AOB')","umax"]
    ammonia = ammonia * mu / (k + ammonia) * (DO / (DO + k_O))
    # logger.debug("ammonia_oxidation:%.2f"%ammonia)
    return ammonia


def nitrite_oxidation(DO, VSS, nitrite):
    mu = NITRITE_OXIDATION_MU * VSS/60
    k = NITRITE_OXIDATION_K
    k_O = NITRITE_OXIDATION_K_O
    # mu_O=kinetics_data.at["('SOUR', 'NOB')","umax"]
    nitrite = nitrite * mu / (k + nitrite) * (DO / (DO + k_O))
    # logger.debug("nitrite_oxidation:%.2f"%nitrite)
    return nitrite


def anammox(VSS, nitrite, inhibition_factor=0.75):
    mu = ANAMMOX_MU * VSS/60
    k = ANAMMOX_K
    inhibition_factor = inhibition_factor
    nitrite = nitrite * mu / (k + nitrite) * inhibition_factor
    # logger.debug("anammox:%.2f"%nitrite)
    return nitrite


def nitrate_denitrification(DO, VSS, nitrate,proportion, conditon='anaerobic'):
    mu = DE_NITRATE_MAX_RATE*VSS

    k = NITRATE_DENITRIFICATION_K
    k_O = 0.2
    if conditon == 'anaerobic':
        nitrate = nitrate * mu / (k + nitrate)
    else:
        nitrate = nitrate * mu / (k + nitrate) * (k_O / (k_O + DO))

    # logger.debug("nitrate_denitrification:%.2f"%nitrate)
    return nitrate


def nitrate_respiration(DO, VSS, nitrate, proportion, conditon='anaerobic'):
    mu = NITRATE_RESPIRATION_MAX_RATE * VSS
    k = NITRATE_RESPIRATION_K
    k_O = 0.2
    if conditon == 'anaerobic':
        nitrate = nitrate * mu / (k + nitrate)
    else:
        nitrate = nitrate * mu / (k + nitrate) * (k_O / (k_O + DO))

    # logger.debug("nitrate_denitrification:%.2f"%nitrate)
    return nitrate


def nitrite_denitrification(DO, VSS, nitrite, conditon='anaerobic'):
    mu = DE_NITRITE_MAX_RATE * VSS
    k = NITRITE_DENITRIFICATION_K
    k_O = 0.2
    if conditon == 'anaerobic':
        nitrite = nitrite * mu / (k + nitrite)
    else:
        nitrite = nitrite * mu / (k + nitrite) * (k_O / (k_O + DO))
    # logger.debug("nitrite_denitrification:%.2f"%nitrite)
    return nitrite


def COD_oxidation(DO, VSS, Islimited=False):
    '''
    COD好氧氧化过程
    :param DO:
    :param Islimited: COD底物是否受限
    :return:
    '''
    mu = COD_OXIDATION_RATE_BY_SOUR * VSS
    k = COD_OXIDATION_K
    if not Islimited:
        COD = DO * mu / (DO + k)

    return COD


def COD_anaerobic_digestion(VSS, COD, Islimited=False):
    mu_constant = COD_ANAEROBIC_DIGESTION_RATE
    if not Islimited:
        COD = mu_constant * VSS
    return COD

def calulate_microaerobic_condition(COD_input, ammonia_input, DO, VSS, COD, ammonia, nitrite, nitrate,
                                    COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR,
                                    ANAMMOX_INHIBITION_FACTOR):
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


        anammox_nitrate = anammox(VSS, nitrite,inhibition_factor=ANAMMOX_INHIBITION_FACTOR) * 0.26 / 1.32
        nitrate_transfer['anammox'] += anammox_nitrate
        NOB_nitrate = nitrite_oxidation(DO, VSS, nitrite)
        nitrate_transfer["nitrite_oxidation"] += NOB_nitrate

        nitrate2n2 = COD_current * COD_FRAC_FOR_DE_NITRATE(DO,VSS,nitrite,nitrate,IsMicroaerobic=True,
                                                           COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR=COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR) / DE_NITRATE_COD
        nitrate2nitrite = COD_current * COD_FRAC_FOR_DE_NITRATE_RESPIRATION(DO,VSS,nitrite,nitrate,IsMicroaerobic=True,
                                                                            COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR=COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR) / NITRATE_RESPIRATION_COD

        DE_nitrate = nitrate2n2 + nitrate2nitrite
        nitrate_transfer["nitrate_denitrification"] += nitrate2n2
        nitrate_transfer["nitrate_respiration"] += nitrate2nitrite
        next_nitrate = nitrate + anammox_nitrate + NOB_nitrate - DE_nitrate - diluent(nitrate)


        AOB_nitrite = ammonia_oxidation(DO, VSS, ammonia)
        anammox_nitrite = anammox(VSS, nitrite,inhibition_factor=ANAMMOX_INHIBITION_FACTOR)
        nitrite_transfer['anammox'] += anammox_nitrite
        NOB_nitrite = nitrite_oxidation(DO, VSS, nitrite)
        nitrite_transfer["nitrite_oxidation"] += NOB_nitrite
        DE_nitrite = COD_current * COD_FRAC_FOR_DE_NITRITE(DO,VSS,nitrite,nitrate,IsMicroaerobic=True,
                                                           COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR=COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR) / DE_NITRITE_COD  # 按照比例使用COD
        nitrite_transfer['nitrite_denitrification'] += DE_nitrite
        nitrite_transfer["nitrate_respiration"] += nitrate2nitrite
        next_nitrite = nitrite + AOB_nitrite + nitrate2nitrite - anammox_nitrite - NOB_nitrite - DE_nitrite - diluent(
            nitrite)
    else:
        # COD充足，不是限制性条件

        next_COD = COD_input + COD - total_COD_next_reaction
        anammox_nitrate = anammox(VSS, nitrite,inhibition_factor=ANAMMOX_INHIBITION_FACTOR) * 0.26 / 1.32
        nitrate_transfer['anammox'] += anammox_nitrate
        NOB_nitrate = nitrite_oxidation(DO, VSS, nitrite)
        nitrate_transfer["nitrite_oxidation"] += NOB_nitrate
        nitrate2n2 = nitrate_rate
        nitrate2nitrite = nitrate_respiration_rate
        DE_nitrate = nitrate2nitrite + nitrate2n2
        nitrate_transfer["nitrate_denitrification"] += nitrate2n2
        nitrate_transfer["nitrate_respiration"] += nitrate2nitrite
        next_nitrate = nitrate + anammox_nitrate + NOB_nitrate - DE_nitrate - diluent(nitrate)


        AOB_nitrite = ammonia_oxidation(DO, VSS, ammonia)
        anammox_nitrite = anammox(VSS, nitrite,inhibition_factor=ANAMMOX_INHIBITION_FACTOR)
        nitrite_transfer['anammox'] += anammox_nitrite
        NOB_nitrite = nitrite_oxidation(DO, VSS, nitrite)
        nitrite_transfer["nitrite_oxidation"] += NOB_nitrite
        DE_nitrite = nitrite_rate
        nitrite_transfer['nitrite_denitrification'] += DE_nitrite
        nitrite_transfer["nitrate_respiration"] += nitrate2nitrite
        next_nitrite = nitrite + AOB_nitrite + nitrate2nitrite - anammox_nitrite - NOB_nitrite - DE_nitrite - diluent(
            nitrite)

    anammox_ammonia = anammox(VSS, nitrite,inhibition_factor=ANAMMOX_INHIBITION_FACTOR) / 1.32
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
             accuracy=0.01, figure=False,
            COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR=1.0,
             ANAMMOX_INHIBITION_FACTOR=0.75,
             ) -> list:
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

    # 用来跟踪底物的浓度随时间的变化，默认只跟踪三氮的变化
    if not TrackList:
        TrackList = ['ammonia', 'nitrite', 'nitrate']
    # 设为带默认值的字典类型数据
    CurrentPeriod_data = defaultdict(list)
    LastPeriod_data = defaultdict(list)
    CurrentPeriod_Mean = defaultdict(float)
    LastPeriod_Mean = defaultdict(float)
    # 结果保存路径
    if not result_path:
        result_path = './result'
    # 当存在result路径，需要提前建立
    if not os.path.exists(result_path):
        try:
            os.makedirs(result_path)
        except FileExistsError:
            pass

    # 过程仿真数据输出到csv文件中,文件的命名规则
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
    logger.debug(f"{ammonia_input}, {COD_input}")
    # logger.debug("ammonia", 'nitrite', 'nitrate')

    # 这里是一个曝气周期的情况。1代表曝气，0代表不曝气
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
                                                                                     COD,ammonia, nitrite, nitrate,
                                                                                     COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR,
                                                                                     ANAMMOX_INHIBITION_FACTOR)
            else:
                for sec in range(0, 60):  # 厌氧状态

                    COD, nitrite, nitrate, ammonia = calculate_anaerobic_condition(COD_input, ammonia_input, 0, VSS,
                                                                                   COD,
                                                                                   ammonia, nitrite, nitrate)

            logger.debug("{0:.2f}  {1:.2f}  {2:.2f}  {3:.2f}".format(ammonia, nitrite, nitrate, COD))

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

        logger.debug(
            "ammonia removal:\nammonia oxidation:{0:.4f}\tanammox:{1:.4f}".format(ammonia_transfer["ammonia_oxidation"],
                                                                                  ammonia_transfer['anammox']))
        logger.debug("nitrite removal:\nnitrite_oxidation:{0:.4f}\tanammox:{1:.4f}\tnitrite_denitrification{2:.4f}".format(
            nitrite_transfer["nitrite_oxidation"], nitrite_transfer['anammox'],
            nitrite_transfer['nitrite_denitrification']))
        logger.debug("nitrate removal:\nnitrate_denitrification{0:.4f}".format(nitrate_transfer['nitrate_denitrification']))
        logger.debug("nitrate generate:\nnitrate_oxidation{0:.2f}\tanammox:{1:.02f}".format(
            nitrate_transfer["nitrite_oxidation"], nitrate_transfer['anammox']))

        for substrate in TrackList:
            CurrentPeriod_Mean[substrate] = np.mean(CurrentPeriod_data[substrate])
        # 判断是否反应稳定
        if not LastPeriod_Mean:
            LastPeriod_Mean = CurrentPeriod_Mean.copy()
        else:


            result_flags = [
                fabs(CurrentPeriod_Mean[substrate] - LastPeriod_Mean[substrate]) <= accuracy for
                substrate in TrackList]
            # 仿真时间不足14小时，且存在一个监控指标误差大于0.001，则再次计算
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
    #logger.debug(f"{params}")
    return simulate(ammonia_init=params["ammonia_init"],
                   nitrite_init=params["nitrite_init"],
                   nitrate_init=params["nitrate_init"],
                   COD_init=params["COD_init"],
                   COD_input=params["COD_input"],
                   ammonia_input=params["ammonia_input"],
                   DO=params["DO"],
                    VSS=params["VSS"],
                    result_path=params["result_path"],
                    figure=params["figure"],
                    COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR=params["COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR"],
                    ANAMMOX_INHIBITION_FACTOR=params["ANAMMOX_INHIBITION_FACTOR"],
                    anaerobic_period=params["anaerobic_period"],
                    microaerobic_period=params["microaerobic_period"]
                    )


def simulate_params(experimental_datas:pd.DataFrame,**kwargs):
    '''

    :param experimental_datas:
    :param kwargs:simulate（）其他的关键字参数，均应为list变量
    :return:
    '''
    data_length=experimental_datas.shape[0]
    # logger.debug(f"{kwargs}")
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

def simulate_coroutines(simulate_proxy, simulate_params_list):
    # 4核运行
    with futures.ProcessPoolExecutor(max_workers=4) as executor:
        for result in executor.map(simulate_proxy, simulate_params_list):

            with open(os.path.join(simulate_params_list[0]["result_path"], 'simulation_data.csv'), mode='a', encoding='utf-8') as f:
                f.write(','.join([str(x) for x in [*result]]) + '\n')



if __name__ == '__main__':
    # 命令行解析
    parser = argparse.ArgumentParser()
    # 结果保存路径
    parser.add_argument('--result_path', '-p',type=str,default="./simulation_period_effect", help='problem class')
    # 是否是追加模式，一般追加模式，是防止程序卡死，再次模拟时重新运行用的。
    parser.add_argument('--append', '-a', type=bool,default=False, help='problem class')
    args = parser.parse_args()
    # 保存路径
    result_path=args.result_path

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
    microaerobic_period_list = [24]
    anaerobic_period_list=[0]
    simulate_list=[(x,y) for x in microaerobic_period_list for y in anaerobic_period_list]

    # 如果之前选择了追加模式，那么就在这里解析最后一条仿真文件夹中包含的仿真参数
    print(args.append)
    if args.append:
        max_time = 0
        for file in os.listdir(result_path):
            current_file_creat_time=os.path.getctime(os.path.join(result_path, file))
            if current_file_creat_time>= max_time:
                max_time=current_file_creat_time
        latest_file_created=file
        # 解析出来用的最后文件的仿真参数,匹配浮点数和整数
        anaerobic_period_str=re.findall('\d+.?\d+',latest_file_created)
        begin_pos=anaerobic_period_list.index(float(anaerobic_period_str))
        # 将旧的已经仿真的内容去除，然后生成新的仿真参数列表
        anaerobic_period_list=anaerobic_period_list[begin_pos:]
        dir_need_del = os.path.join(result_path, latest_file_created)
        for file in os.listdir(dir_need_del):
            os.remove(os.path.join(dir_need_del, file))
        # 最后删除空的文件夹
        os.rmdir(dir_need_del)

    #变成tqdm对象，方便输出进度条
    simulate_list = tqdm(simulate_list)
    for microaerobic_period,anaerobic_period in simulate_list:
        # 设定进度条的格式
        simulate_list.set_description(f"Processing:microaerobic_period:{microaerobic_period} and anaerobic_period:{anaerobic_period}")
        # 仿真文件，储存的目录:
        simulate_root_path = os.path.join(result_path,f"microaerobic_period_{microaerobic_period}&anaerobic_period_{anaerobic_period}")
        simulate_root_path_list = [simulate_root_path] * simulate_datas.shape[0]
        microaerobic_period = [microaerobic_period] * simulate_datas.shape[0]
        anaerobic_period=[anaerobic_period]*simulate_datas.shape[0]
        cod_anaerobic_digestion_inhibition_factor=0
        anammox_inhibition_factor=0.75
        COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR=[cod_anaerobic_digestion_inhibition_factor]*simulate_datas.shape[0]

        ANAMMOX_INHIBITION_FACTOR=[anammox_inhibition_factor]*simulate_datas.shape[0]
        simulate_params_list = simulate_params(simulate_datas, DO=DO, VSS=VSS_datas,
                                               result_path=simulate_root_path_list,
                                               COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR=COD_ANAEROBIC_DIGESTION_INHIBITION_FACTOR,
                                               ANAMMOX_INHIBITION_FACTOR=ANAMMOX_INHIBITION_FACTOR,
                                               anaerobic_period=anaerobic_period,
                                               microaerobic_period=microaerobic_period,
                                               )
        # 4核运行

        simulate_coroutines(simulate_proxy, simulate_params_list)