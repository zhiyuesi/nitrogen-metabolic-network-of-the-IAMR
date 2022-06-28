# -*- coding: utf-8 -*-
import numpy as np
import geatpy as ea # import geatpy
import importlib
import argparse

import pandas as pd
import os
import matplotlib.pyplot as plt

def run_GA(problem,result_path):
    '''
    运行GA脚本，求解相应问题的动力学参数
    :param problem:
    :param result_path:
    :return:
    '''

    """===============================实例化问题对象==========================="""
    #ea.mutbin()
    problem = problem()  # 生成问题对象
    """=================================种群设置==============================="""
    Encoding = 'BG'  # 编码方式
    NIND = 100  # 种群规模
    Field = ea.crtfld(Encoding, problem.varTypes, problem.ranges, problem.borders)  # 创建区域描述器，自定义的精度为10^-4
    population = ea.Population(Encoding, Field, NIND)  # 实例化种群对象（此时种群还没被初始化，仅仅是完成种群对象的实例化）
    """===============================算法参数设置============================="""
    # 增强精英保留的遗传算法模板
    myAlgorithm = ea.soea_SEGA_templet(problem, population)  # 实例化一个算法模板对象
    myAlgorithm.MAXGEN = 100  # 最大进化代数
    myAlgorithm.drawing = 0  # 不进行绘图

    """==========================调用算法模板进行种群进化======================="""
    # 写入每个问题每一次计算的动力学最优信息
    current_save_path = result_path
    if not os.path.exists(current_save_path):
        os.mkdir(current_save_path)
    with open(os.path.join(current_save_path,'best_var.csv'), mode='w', encoding='utf-8') as f:
        f.write('μmax,Ks,best_ObjV\n')
    # 保存程序的运行目录的位置
    run_path = os.getcwd()
    for i in range(10):  # 每一个运行十次
        [population, obj_trace, var_trace] = myAlgorithm.run()  # 执行算法模板

        os.chdir(current_save_path)
        population.save()  # 把最后一代种群的信息保存自定义的路径中去
        plt.close()
        # 输出结果
        trace = pd.DataFrame(obj_trace, columns=['mean_ObjV', 'Best_ObjV'])
        trace.to_csv('Result/ObjV.csv', encoding='utf-8', index=False)
        os.rename('Result', 'Result' + str(i))
        best_gen = np.argmin(problem.maxormins * obj_trace[:, 1])  # 记录最优种群个体是在哪一代
        # mean_ObjV =  obj_trace[:, 0] # 记录每一次的平均目标函数值

        best_ObjV = obj_trace[best_gen, 1]
        print('最优的目标函数值为：%s' % (best_ObjV))
        print('最优的控制变量值为：')
        with open('best_var.csv', 'a', encoding='utf-8') as f:

            for i in range(var_trace.shape[1]):
                print(var_trace[best_gen, i])
                f.write(str(var_trace[best_gen, i]) + ',')
            f.write(str(best_ObjV))
            f.write('\n')
        print('有效进化代数：%s' % (obj_trace.shape[0]))
        print('最优的一代是第 %s 代' % (best_gen + 1))
        print('评价次数：%s' % (myAlgorithm.evalsNum))
        print('时间已过 %s 秒' % (myAlgorithm.passTime))
        os.chdir(run_path)#改回到原来的运行程序所在的路径

    with open("Optimization.csv",mode="a+") as f:
        f.write(current_save_path+"\t")

        # 查找每个问题类的最优解
        cal_result = pd.read_csv(os.path.join(current_save_path, "best_var.csv"), encoding="utf-8", header=0, index_col=None)
        opt_obj=cal_result.loc[cal_result["best_ObjV"].idxmin()]
        # 将其转换为字符串直接进行拼接
        f.write("\t".join(map(str,opt_obj)))
        f.write("\n")
    #with open(os.path.join(problem,""))
if __name__ == '__main__':
    """===============================配置命令行参数==========================="""
    parser = argparse.ArgumentParser()
    # 动态添加问题类
    parser.add_argument('--problems', '-p', nargs="+",type=str, help='problem class')
    args = parser.parse_args()
    print(args)
    for problem in args.problems:
        problem_class=importlib.import_module(problem).MyProblem
        run_GA(problem_class,problem)

