# -*- coding: utf-8 -*-
'''
@Author  :  月司
@Email   :  815603884@qq.com
@File    :  确定anammox的最佳氧半抑制浓度.py
@version :  1.0
@Desc    : 
'''

lolims, uplims=0.65,0.75

def inhibiton_eq(k_o):
    DO=0.2
    return k_o/(k_o+DO)
epsilon=1e-3
for k_o in range(1500):
    k_o=k_o/1000
    f=inhibiton_eq(k_o)
    if abs(lolims-f)<=epsilon:
        break
k_o_lolims=k_o
for k_o in range(1500):
    k_o=k_o/1000
    f=inhibiton_eq(k_o)
    if abs(uplims-f)<=epsilon:
        break

k_o_uplims=k_o
print(f"anammox的DO半抑制常数在{k_o_lolims}-{k_o_uplims}之间")
