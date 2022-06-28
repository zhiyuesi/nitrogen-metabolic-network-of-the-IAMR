# -*- coding: utf-8 -*-
'''

@Author  :  月司
@Email   :  815603884@qq.com
@File    :  添加仿真的最大值和最小值信息.py
@version :  1.0
@Desc    : 
'''
import pandas as pd
from collections import namedtuple
import os

# 绘制拟合的二维图
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

simulation_root='./result/'
file_list=[]

File = namedtuple("_file","filename,created")

for file in os.listdir(simulation_root):


    if not file=="simulation_data.csv" and os.path.splitext(file)[1]==".csv":

        current_file_creat_time = os.path.getctime(os.path.join(simulation_root, file))
        file_list.append(File(file,current_file_creat_time))
#
file_list=file_list[-14:]
file_list=sorted(file_list,key=lambda file:file.created)

simulation_data = pd.read_csv('./result/simulation_data.csv', encoding='utf-8', header=None)
simulation_data.columns=["ammonia_input","cod_input","ammonia_mean","nitrite_mean","nitrate_mean"]
for column in [x+"_"+y for x in ["ammonia","nitrite","nitrate"] for y in ["min","max"] ]:
    simulation_data[column]=None
substance_list=["ammonia","nitrite","nitrate"]
for index,file in enumerate(file_list):
    substance_trace_data=pd.read_csv(os.path.join(simulation_root,file.filename),header=0,index_col=0)
    substance_trace_data=substance_trace_data[-24:]
    for substance in substance_list:
        # 找到每种物质的最大值和最小值
        simulation_data.at[index,substance+"_min"]=substance_trace_data[substance].min()
        simulation_data.at[index,substance+"_max"]=substance_trace_data[substance].max()

simulation_data.to_csv(os.path.join(simulation_root,"simulation_data_mean_min_max.csv"),index=False,encoding="utf-8")
