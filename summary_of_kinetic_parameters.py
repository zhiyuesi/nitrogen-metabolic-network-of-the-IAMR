# -*- coding: utf-8 -*-
'''
@Author  :  月司
@Email   :  815603884@qq.com
@File    :  summary_of_kinetic_parameters.py
@version :  1.0
@Desc    : 
'''

import os
import pandas as pd

# 读取生物量数据
VSS = pd.read_csv('VSS.csv', header=0, index_col=[0, 1, 2])

path_root = ["Nitrogen_Substance", "SOUR"]
data_result = []
SOUR_control = []
for path in path_root:
    for child_path in os.listdir(path):
        current_path = os.path.join(path, child_path)
        if os.path.isdir(current_path):
            try:
                current_result_data = pd.read_csv(os.path.join(current_path, "Optimization.csv"), delimiter="\t",
                                                  header=None, index_col=0, encoding='UTF-8')
                for problem in current_result_data.index:
                    assay_ID = problem.partition("_")[-1]  # 去掉Myproblem
                    assay_ID = assay_ID.rpartition("_")[0]  # 去掉re

                    if assay_ID:
                        data_result.append([path, child_path, assay_ID,
                                            current_result_data.at[problem, 1],
                                            current_result_data.at[problem, 2],
                                            VSS.at[(path, child_path, assay_ID), "VSS"]])

            except:
                if os.path.exists(os.path.join(current_path, "control.csv")):
                    control_data = pd.read_csv(os.path.join(current_path, "control.csv"), delimiter="\t",
                                               header=None, index_col=0)
                    for assay_ID in control_data.index:
                        SOUR_control.append([path, child_path, assay_ID,
                                             control_data.at[assay_ID, 1],
                                             VSS.at[(path, child_path, assay_ID), "VSS"]])

data_result = pd.DataFrame(data_result, columns=["type", "assay_target", "assay_id", "umax", "ks", "VSS"])
SOUR_control = pd.DataFrame(SOUR_control, columns=["type", "assay_target", "assay_id", "umax", "VSS"])

# 变成multi-index
data_result.set_index(["type", "assay_target", "assay_id"], inplace=True)
SOUR_control.set_index(["type", "assay_target", "assay_id"], inplace=True)
data_result.to_json("assay_kinetics_details.json", orient="index",double_precision=3)
SOUR_control.to_json("SOUR_control_details.json", orient="index",double_precision=3)
# 计算平均生物量的动力学参数
data_result["umax"]=data_result["umax"]/data_result["VSS"]

# sour空白试验的平均速率。
SOUR_control["umax"]=SOUR_control["umax"]/SOUR_control["VSS"]
SOUR_control_rate=SOUR_control["umax"].mean()
SOUR_data=data_result.loc["SOUR"]
SOUR_data["umax"]=SOUR_data["umax"]-SOUR_control_rate
# 寻找每个实验的动力学参数。
mean_result = data_result.groupby(level=[0, 1]).mean()
# 这个都是单位VSS的数据了，所以需要把VSS列去掉。
mean_result.drop(columns=["VSS"],inplace=True)
mean_result.to_json("kinetics.json", orient="index",double_precision=3)
