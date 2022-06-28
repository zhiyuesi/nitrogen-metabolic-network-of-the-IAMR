# -*- coding: utf-8 -*-
'''
@Author  :  月司
@Email   :  815603884@qq.com
@File    :  calculate_nitrogen_conversion.py
@version :  1.0

@Desc    :  计算每个脱氮路径的脱氮贡献
'''

import os
import json
from collections import namedtuple

result_path="./result"

# 因为保存的json开头是进水的COD和NH4+-N，所以文件的默认排序跟时间顺序是不同的
order_filelist=[]
for file in os.listdir(result_path):
    if file.rpartition(".")[-1]=="json":
        order_filelist.append(file)

# 进行排序
def sort_by_time(item:str):
    full_path=os.path.join(result_path,item)
    mtime = os.stat(full_path).st_mtime
    return  mtime

order_filelist.sort(key=sort_by_time)

result_list=[]
for file in order_filelist:
    # 这个rsult_dict保存不同底物在不同路径上的转移量。
    result_dict={}
    #注意，这些json文件的转移量都是正值
    result=json.loads(open(os.path.join(result_path,file)).read())
    #有关ammonia的转化过程
    result_dict["ammonia"]=result[0]
    result_dict["nitrite"] = result[1]
    result_dict["nitrate"] =result[2]
    result_list.append(result_dict)



# 计算anammox和反硝化的脱氮贡献
import numpy as np
anammox_list=[]
nitrite_list=[]
nitrate_list=[]
de_list=[]
for item in result_list:
    anammox_contribution=item["ammonia"]["anammox"]+item["nitrite"]["anammox"]-item["nitrate"]["anammox"]
    #亚硝酸盐反硝化
    nitrite_contribution=item["nitrite"]["nitrite_denitrification"]
    #硝酸盐反硝化
    nitrate_contribution=item["nitrate"]["nitrate_denitrification"]
    # 硝酸盐呼吸
    nitrate_respiration_contribution = item["nitrate"]["nitrate_respiration"]
    # 这里计算的总共去除的N量。
    total_removed=anammox_contribution+nitrate_contribution+nitrite_contribution
    # anammox去除的氮量与总去除量的比值。
    anammox_list.append(anammox_contribution/total_removed)
    nitrite_list.append(nitrite_contribution/total_removed)
    nitrate_list.append(nitrate_contribution/total_removed)

print(f"anammox的总的脱氮效率{np.mean(anammox_list)}")
print(f"亚硝酸盐反硝化的总的脱氮效率{np.mean(nitrite_list)}")
print(f"硝酸盐反硝化的总的脱氮效率{np.mean(nitrate_list)}")

# anammox的亚硝酸盐利用率：就是氧化得到的亚硝酸盐有哪些是通过anammox的途径去除的
anammox_nitrite_list=[]
for item in result_list:
    anammox_contribution=item["nitrite"]["anammox"]

    total_nitrite=item["ammonia"]["ammonia_oxidation"]
    anammox_nitrite_list.append(anammox_contribution/total_nitrite)

print(f"anammox的亚硝酸盐利用率:{np.mean(anammox_nitrite_list)}")

# 硝酸盐积累的量
delta_nitrate=[]
for item in result_list:
    nitrate_rest=item["nitrate"]["anammox"]+item["nitrate"]["nitrite_oxidation"]\
                 -item["nitrate"]["nitrate_denitrification"]-item["nitrate"]["nitrate_respiration"]
    total_transferd=item["ammonia"]["anammox"]+item["ammonia"]["ammonia_oxidation"]

    delta_nitrate.append(nitrate_rest/total_transferd)

print(f"硝酸盐积累量:{np.mean(delta_nitrate)}")

# anammox生成硝酸盐所占的百分比
nitrate_anammox_list=[]
for item in result_list:
    nitrate_anammox=item["nitrate"]["anammox"]
    total_transferd=item["nitrate"]["anammox"]+item["nitrate"]["nitrite_oxidation"]
    nitrate_anammox_list.append(nitrate_anammox/total_transferd)
print(f"anammox生成的硝酸盐的百分比:{np.mean(nitrate_anammox_list)}")

# 生成一个关于去除路径的平均值。
import copy

# 用实际的NH4+-N输入标准化
for index,item in enumerate(result_list):
    filename=order_filelist[index]
    base=float(filename.split("_")[1])/8*(24/60)
    # item保存的是一天的氮转移速率的字典，key为NH4+-N等物质，value是包含其所有转移途径转移量的字典。
    for key, value in item.items():
        if isinstance(value, dict):
        # key2是具体的转移途径，value2是转移量
            for key2,value2 in value.items():
                value[key2]=value2/base
        else:
            item[key]=value/base
print(f"各天标准化的数据：{result_list}")

# 保存平均去除路径的字典
transfer_dict=copy.deepcopy(result_dict)
TRANSFER=namedtuple("_transfer",["mean","min","max"])
for key,value in transfer_dict.items():
    for key, value in item.items():
        if isinstance(value, dict):
            for key2,value2 in value.items():
                temp_list=[]
                for day_data in result_list:

                    temp_list.append(day_data[key][key2])
                transfer_dict[key][key2]=TRANSFER(np.mean(temp_list),np.min(temp_list),np.max(temp_list))
        else:
            temp_list = []
            for day_data in result_list:
                temp_list.append(day_data[key])
            transfer_dict[key][key2]=TRANSFER(np.mean(temp_list),np.min(temp_list),np.max(temp_list))
print(f"平均数据{transfer_dict}")


