#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author  : 月司
# @Email   : 815603884@qq.com
# @File    : demo_C_s_AOB.py
#@version  : 1.0

import matplotlib.pyplot as plt
import pandas as pd

#latex 和times new roman
plt.rcParams['font.sans-serif']=['times']
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family']=['Times New Roman']


# 载入动力学文件
kinetics_data=pd.read_json("../kinetics.json",orient="index")

# 载入生物量文件
VSS = pd.read_csv('../VSS.csv', header=0, index_col=[0, 1, 2])
#
VSS_list=[
	VSS.at[("Composite_Processes","ANAMMOX&Denitrification_NO2","A&De_nitrite1"),"VSS"],
	VSS.at[("Composite_Processes","ANAMMOX&Denitrification_NO2","A&De_nitrite2"),"VSS"],
	VSS.at[("Composite_Processes","ANAMMOX&Denitrification_NO2","A&De_nitrite3"),"VSS"],
		  ]


def anammox(nitrite,sample_VSS,inhibition_factor=0.75):
	mu = kinetics_data.at["('Nitrogen_Substance', 'ANAMMOX')","umax"]*sample_VSS
	k = kinetics_data.at["('Nitrogen_Substance', 'ANAMMOX')","ks"]
	inhibition_factor=inhibition_factor
	nitrite=nitrite * mu / (k + nitrite)*inhibition_factor
	#print("anammox:%.2f"%nitrite)
	return nitrite

def nitrite_denitrification(nitrite,sample_VSS):
	mu = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO2')","umax"]*sample_VSS
	k = kinetics_data.at["('Nitrogen_Substance', 'Denitrification_NO2')","ks"]
	nitrite=nitrite * mu / (k + nitrite)
	#print("nitrite_denitrification:%.2f"%nitrite)
	return nitrite

if __name__ == '__main__':
	fig_width = 18
	fig_height = 6
	width = fig_width / 2.54
	height = fig_height / 2.54
	#dpi = 600
	dpi=96
	fig,ax=plt.subplots(1,3,figsize=(width,height),dpi=dpi,sharey=True)
	fig.subplots_adjust(wspace=0)
	fig.subplots_adjust(0.2, 0.2, 0.95, 0.85)

	data_1=pd.read_csv("./data/A&De_nitrite1.csv",header=0,index_col=None)
	data_2 = pd.read_csv("./data/A&De_nitrite2.csv", header=0, index_col=None)
	data_3 = pd.read_csv("./data/A&De_nitrite3.csv", header=0, index_col=None)

	# 实验数据
	sample_list=[data_1,data_2,data_3]
	#
	# sample是一个dataframe的实验数据
	for idx,sample in enumerate(sample_list):
		# 去掉第一个数据，因为可能由于吸附等因素，导致不准。
		C_s_list= list(sample.loc[1:,"NO2--N"])
		delta_time = []
		# 因为去掉第一个数据，所以计算时间间隔也要从1开始计算
		for index in sample.index[2:]:
			delta_time.append(sample.at[index, "Mins"] - sample.at[index - 1, "Mins"])
		print(delta_time)

		cal_data = []

		delta = 0
		# C_s_list[0]是实验数据的第二个点，从第二个点一直模拟到实验结束。
		C_s=C_s_list[0]
		for delta in delta_time:
			for i in range(delta):
			#这里是每分钟模拟一次
				C_s = C_s - anammox(C_s,VSS_list[idx],inhibition_factor=1.0)-nitrite_denitrification(C_s,VSS_list[idx])
			cal_data.append(C_s)
		#因为没有办法预测首个点，所以将初始点插入到预测列表的开头，让其长度等于实际数据数据
		cal_data.insert(0,C_s_list[0])
		print(C_s_list)
		print(cal_data)
		ax[idx].scatter(sample.loc[1:,'Mins'],C_s_list,c='red')
		ax[idx].plot(sample.loc[1:,'Mins'],cal_data)
		ax[idx].set_xlabel(r"Time (min)", ha="center",
			   va="center")
	ax[0].set_ylabel(r"$\rm NO_2^--N$ (mg/L)", ha="center",
			   va="center")
	plt.savefig('Competition_for_nitrite.jpg')
	plt.show()
