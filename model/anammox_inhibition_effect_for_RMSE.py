# -*- coding: utf-8 -*-
'''
@Author  :  月司
@Email   :  815603884@qq.com
@File    :  初步探究anammox最佳抑制比.py
@version :  1.0
@Desc    : 
'''

import pandas as pd
import os
from collections import defaultdict
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



RMSE_list = defaultdict(list)
for factor in [x/100 for x in range(0,105,5)]:
    # 仿真文件储存的根目录:
    simulate_root_path = f"./inhibition/anammox{factor}&COD{0.0}"
    print(f"抑制因子为{factor}时，",end="\t")

    simulation_data = pd.read_csv(os.path.join(simulate_root_path, 'simulation_data.csv'), encoding='utf-8',
                                  header=None)
    simulation_data.columns = ["COD_Input", "Ammonia_Input", "Ammonia_OUT", "Nitrite_OUT", "Nitrate_OUT"]
    show_list = ["Ammonia_OUT", "Nitrite_OUT", "Nitrate_OUT"]

    for (index, substrate) in enumerate(show_list):
        # F代表真实的数据
        experimental_data_label = substrate.split("_")[0] + "_F"
        substrate_experimental_data=experimental_data[experimental_data_label][-14:]
        substrate_simulation_data=list(simulation_data[substrate])[-14:]
        RMSE=0
        for experiment,simulate in zip(substrate_experimental_data,substrate_simulation_data):
            # 亚硝酸盐，因为仿真数据总是大于真实数据。所以反过来求解
            if substrate == "Nitrite_OUT":
                RMSE += ((simulate-experiment) / simulate) ** 2
            else:
                RMSE += ((experiment - simulate) / experiment) ** 2
        RMSE=RMSE/(len(substrate_experimental_data)-1)
        RMSE_list[substrate].append(RMSE)
        print(f"{substrate}的平均值为{round(np.mean(substrate_simulation_data),2)}，的RMSE为{round(RMSE,6)}",end="\t")
    print("\n")

fig, ax = plt.subplots(3, 1, sharex=True)
fig.subplots_adjust(hspace=0)

show_list = ["Ammonia_OUT", "Nitrite_OUT", "Nitrate_OUT"]
for (index, substrate) in enumerate(show_list):
    ax[index].plot([x/100 for x in range(0,105,5)],RMSE_list[substrate], color='blue', marker='o', label=substrate)
    ax[index].legend()
plt.legend()
plt.show()