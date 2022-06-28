# -*- coding: utf-8 -*-
'''
@Author  :  月司
@Email   :  815603884@qq.com
@File    :  summary_inhibit_factors_RMSE.py
@version :  1.0
@Desc    : 
'''


import pandas as pd
import os
from collections import  defaultdict
import matplotlib.pyplot as plt
import  numpy as np

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

# 用于仿真的实验数据
simulate_datas = experimental_data.loc[42:,
                 ['Ammonia_OUT', 'Nitrite_OUT', 'Nitrate_OUT', 'COD_OUT', 'COD_removal', 'TKN']]
simulate_root_path=f"./inhibition"
# anammox抑制因子列表,0-1之间均匀的21个点
anammox_inhibition_factor_list = [x / 100 for x in range(0, 105, 5)]
# COD抑制因子，0-0.1十-点，0.1-1是9个点。
cod_anaerobic_digestion_inhibition_factor_list = [x / 100 for x in range(10)] + [x / 10 for x in range(1, 11)]
combined_inhibition_factor_list = [(x, y) for x in anammox_inhibition_factor_list for y in
                                   cod_anaerobic_digestion_inhibition_factor_list]

# key为底物的名称
RMSE_list = defaultdict(list)
for anammox_inhibition_factor, cod_anaerobic_digestion_inhibition_factor in combined_inhibition_factor_list:

    # 仿真文件储存的根目录
    current_simulate_path = os.path.join(simulate_root_path,f"anammox{anammox_inhibition_factor}&COD{cod_anaerobic_digestion_inhibition_factor}")

# 模型仿真得到的数据
    simulation_data = pd.read_csv(os.path.join(current_simulate_path, 'simulation_data.csv'), encoding='utf-8',
                                  header=None)
    simulation_data.columns = ["COD_Input", "Ammonia_Input", "Ammonia_OUT", "Nitrite_OUT", "Nitrate_OUT"]
    # 这些标签是simulate数据的
    show_list = ["Ammonia_OUT", "Nitrite_OUT", "Nitrate_OUT"]

    for (index, substrate) in enumerate(show_list):
        # F代表真实的数据
        experimental_data_label = substrate.split("_")[0] + "_F"
        substrate_experimental_data=experimental_data[experimental_data_label][-14:]
        substrate_simulation_data=list(simulation_data[substrate])[-14:]
        RMSE=0
        for experiment,simulate in zip(substrate_experimental_data,substrate_simulation_data):
            # 亚硝酸盐，因为仿真数据总是大于真实数据。所以就反过来
            if substrate=="Nitrite_OUT":
                RMSE += ((-experiment + simulate) / simulate) ** 2
            else:
                RMSE+=((experiment-simulate)/experiment)**2
        RMSE=RMSE/(len(substrate_experimental_data)-1)
        RMSE_list[substrate].append([anammox_inhibition_factor,cod_anaerobic_digestion_inhibition_factor,RMSE])

RMSE_data_dict={}
for key in RMSE_list:
    RMSE_data_dict[key]=pd.DataFrame(RMSE_list[key],columns=[
        "anammox_inhibition_factor",
        "cod_anaerobic_digestion_inhibition_factor",
        key
    ])
# 合并DataFrame为一个
key_list=list(RMSE_data_dict.keys())
left_df=RMSE_data_dict[key_list[0]]
for pos in range(1,len(key_list)):
    left_df=pd.merge(left_df,RMSE_data_dict[key_list[pos]])
#保存为到本地文件中
left_df.to_csv("inhibit_factors_RMSE.csv",encoding="utf-8")
