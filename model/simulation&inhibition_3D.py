# -*- coding: utf-8 -*-
'''
@Author  :  月司
@Email   :  815603884@qq.com
@File    :  绘制三维图.py
@version :  1.0
@Desc    : 
'''

import pandas as pd
import os
from collections import  defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import IndexLocator,MultipleLocator,AutoLocator,LinearLocator
from mpl_toolkits.mplot3d import Axes3D
import  numpy as np
from matplotlib.pyplot import  Line2D
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = ['serif']
plt.rcParams['font.serif'] = ['Times New Roman',"Times"]


fig = plt.figure(figsize=(6.0, 6.01),dpi=100)

num_rows = 4
num_cols = 5
# 这里的行和空白的高度
row_height = 2

#空白的高度
# 行空白
row_space_height = 0
# 列空白
col_space_height = 1

# 网络规模,行2*5+2=11，列=2*5=13
grid = (row_height*num_rows + row_space_height, num_cols*row_height+col_space_height)


# 绘制三维图
RMSE_data=pd.read_csv("inhibit_factors_RMSE.csv",header=0,index_col=None,encoding="utf-8")
anammox_inhibition_factor_list=list(RMSE_data["anammox_inhibition_factor"].unique())
cod_anaerobic_digestion_inhibition_factor_list=list(RMSE_data["cod_anaerobic_digestion_inhibition_factor"].unique())
# 生成网格点坐标系统
X,Y = np.meshgrid(anammox_inhibition_factor_list,cod_anaerobic_digestion_inhibition_factor_list)
RMSE_data.set_index(["anammox_inhibition_factor","cod_anaerobic_digestion_inhibition_factor"],inplace=True)
show_list = ["Ammonia_OUT", "Nitrite_OUT", "Nitrate_OUT"]

plot_list=[
    {"loc":(0,num_cols+col_space_height),
     "rowspan":num_rows,
    "colspan":num_cols,},
    {"loc": (num_rows+row_space_height,0),
     "rowspan": num_rows,
     "colspan": num_cols, },
    {"loc": (num_rows+row_space_height,num_cols+col_space_height),
     "rowspan": num_rows,
     "colspan": num_cols, },

]



anammox_inhibition_dir=""
for index,substance in enumerate(show_list):
    #
    Z=np.zeros_like(X)
    for anammox_inhibition_factor, cod_anaerobic_digestion_inhibition_factor in RMSE_data.index:
        X_pos=anammox_inhibition_factor_list.index(anammox_inhibition_factor)
        Y_pos=cod_anaerobic_digestion_inhibition_factor_list.index(cod_anaerobic_digestion_inhibition_factor)
        #要注意，矩阵的行为X,矩阵的列为Y
        Z[Y_pos,X_pos]=RMSE_data.at[(anammox_inhibition_factor,cod_anaerobic_digestion_inhibition_factor),substance]


    ax=plt.subplot2grid(grid, plot_list[index]["loc"],
                     rowspan=plot_list[index]["rowspan"],
                     colspan=plot_list[index]["colspan"],
                        fig=fig,
                     projection='3d')

    surf = ax.plot_surface(X, Y, Z, rstride=1,
                           cstride=1,
                           cmap='coolwarm',
                           # 能提高显示的效果，抗锯齿
                           antialiased=False,
                           shade=False,

                           )


    ax.set_xlabel("AIF", fontsize=9,labelpad=-5 )
    ax.set_ylabel("CIF", fontsize=9,labelpad=-4 )
    ax.set_zlabel("ReMSE", fontsize=9,labelpad=-6)

    ax.tick_params(direction="in",reset=True,pad=-2,length=10)


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

simulation_data = pd.read_csv('./result/simulation_data_mean_min_max.csv', encoding='utf-8', header=0,index_col=False)
simulation_data=simulation_data[-14:]
show_list = ["Ammonia_OUT", "Nitrite_OUT", "Nitrate_OUT"]


simulation_ylim=[(15,55),(0,1.6),(10,48)]
# 给生成刻度生成器的base用的(MultipleLocator)
num_rows = 6
num_cols = 6
# 这里的行和空白的高度，的单位是同一个。
row_height = 2
#空白的高度,直观的解释就是这个空白占画框高度的比例为多少
# 行空白
row_space_height = 2
# 列空白
col_space_height = 2

grid = (row_height*num_rows + row_space_height, num_cols*row_height+col_space_height)
simulation_plot=[
    {"loc":(1,col_space_height),
     "rowspan":row_height,
    "colspan":num_cols-1,},
    {"loc": (row_height+1,col_space_height),
     "rowspan": row_height,
     "colspan": num_cols-1, },
    {"loc": (row_height*2+1,col_space_height),
     "rowspan": row_height,
     "colspan": num_cols-1, },

]

substrate_label= {
"Ammonia":r"$\rm NH_4^+-N$",
"Nitrite":r"$\rm NO_2^--N$",
"Nitrate":r"$\rm NO_3^--N$",
}


simulation_ylabel_locator=[10,0.5,10]
for (index, substance) in enumerate(show_list):
    # F代表真实的数据
    markersize = 4
    lw=1
    ax=plt.subplot2grid(grid, simulation_plot[index]["loc"],
                     rowspan=simulation_plot[index]["rowspan"],
                     colspan=simulation_plot[index]["colspan"],
                    fig=fig
                     )

    ax.tick_params(direction="in")
    plt.subplots_adjust(hspace=0)
    substance=substance.split("_")[0]
    experimental_data_label = substance + "_F"
    ax.plot(range(44, 58),list(experimental_data[experimental_data_label])[-14:], marker='o',markersize=markersize,lw=lw, color='red',
                   label=experimental_data_label)
    ylowerr= simulation_data[substance.lower() + "_mean"] - simulation_data[substance.lower() + "_min"]
    yuperr= simulation_data[substance.lower() + "_max"] - simulation_data[substance.lower() + "_mean"]
    yerr=[ylowerr,yuperr]

    ax.plot(range(44, 58),list(simulation_data[substance.lower() + "_mean"]),
            color='blue', marker='o',ms=markersize,
            label=substance.lower())
    ax.set_ylabel(substrate_label[substance] + "\n(mg/L)", va="bottom", fontdict={
        "fontsize":9
    })
    ax.set_xlim(43.5, 58)
    ax.set_xlabel("Time (day)")
    ax.xaxis.set_major_locator(MultipleLocator(4))
    ylim_min,ylim_max=simulation_ylim[index]
    ax.set_ylim(ylim_min,ylim_max)
    ax.yaxis.set_major_locator(MultipleLocator(simulation_ylabel_locator[index]))
    plt.subplots_adjust(bottom=.0, top=1.0,left=0,right=0.95)
    # ax.legend()


plt.annotate(
    text="(A)",xy=((num_cols)/grid[1],(grid[0]-row_space_height+0.5)/grid[0]),
    xycoords="figure fraction"
)
plt.annotate(
    text="(B)",xy=((num_cols*row_height+1)/grid[1],(grid[0]-row_space_height+0.5)/grid[0]),
    xycoords="figure fraction"
)

plt.annotate(
    text="(C)",xy=((num_cols)/grid[1],(num_rows-0.5)/grid[0]),
    xycoords="figure fraction"
)
plt.annotate(
    text="(D)",xy=((num_cols*row_height+1)/grid[1],(num_rows-0.5)/grid[0]),
    xycoords="figure fraction"
)
plt.subplots_adjust(wspace=.0,hspace=.0)


legend_elements = [
                   Line2D([0], [0], marker='o',lw=lw, color='blue', label='measured value',
                          markerfacecolor='b', markersize=markersize),
                    Line2D([0], [0], marker='o',color='red', lw=lw,markersize=markersize, label='simulated value',),
]

fig.legend(handles=legend_elements,
           frameon=False,
           ncol=2,
            columnspacing = 0.5,
           labelspacing = 0.1,
           loc='upper left',bbox_to_anchor=(0.02,0.98)
           )
plt.show()
plt.savefig("simulation&inhibition_3D.jpg",dpi=600)