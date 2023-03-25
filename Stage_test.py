from API.Stage import Stage
from API.TCPv2 import ThPoint
from API.Stage import GridProfile

from pprint import pprint
import numpy as np
import pandas as pd


import time

#aaa = GridProfile(M=0.414,type=1,in_angle=27.591,out_angle=17.611)

# print(int(301.111))
#print(aaa.get_grids_list())
# df = pd.DataFrame(columns=[
#     'Name',
#     'b',
#     'B',
#     'a',
#     'rel_t',
#     'ins_angle',
#     'Delta_in',
#     'Delta_out',
#     'center_X',
#     'center_Y',
#     'F',
#     'Ixx',
#     'Iyy',
#     'Wxx_back',
#     'Wxx_edge',
#     'Wyy_in_edge',
#     'Wyy_out_edge'
# ])

# df = df.append({
#     'Name':'C-04-19',
#     'b':1,
#     'B':1,
#     'a':1,
#     'rel_t':1,
#     'ins_angle':1,
#     'Delta_in':1,
#     'Delta_out':1,
#     'center_X':1,
#     'center_Y':1,
#     'F':1,
#     'Ixx':1,
#     'Iyy':1,
#     'Wxx_back':1,
#     'Wxx_edge':1,
#     'Wyy_in_edge':1,
#     'Wyy_out_edge':1
# },ignore_index=True)

#df.to_csv('aue.csv')

# df1 = pd.read_excel('ffff.xlsx')
# pprint(df1)

# df2 = df1
# df2['b'] = df2['b']/1000
# df2['B'] = df2['B']/1000
# df2['a'] = df2['a']/1000
# df2['Delta_in'] = df2['Delta_in']/1000
# df2['Delta_out'] = df2['Delta_out']/1000
# df2['center_X'] = df2['center_X']/1000
# df2['center_Y'] = df2['center_Y']/1000

# df2['F'] = df2['center_X']/100/100

# df2['Ixx'] = df2['Ixx']/100/100/100/100
# df2['Iyy'] = df2['Iyy']/100/100/100/100

# df2['Wxx_back'] = df2['Wxx_back']/100/100/100
# df2['Wxx_edge'] = df2['Wxx_edge']/100/100/100
# df2['Wyy_in_edge'] = df2['Wyy_in_edge']/100/100/100
# df2['Wyy_out_edge'] = df2['Wyy_out_edge']/100/100/100

# df2.to_excel('API/grids.xlsx')

# def fuck(dt:dict):
#     for i in dt:
#         print(i)

# fuck({
#     "bla":1,
#     "dra":13,
# })
    
# print(len({
    
# }))
# f = Stage(
#     p_0=1.103,#1.103,
#     t_0=410.7,
#     G_0=185.7,
#     d_hub=2.275-0.133,
#     n = 50,
#     reaction=0.337,
#     alpha_0 = 90,
#     alpha_1eef=10,
#     H_0 = 184.5,#184.5,
#     c_0 = 0,
#     Delta_pr=0.013,
#     kappa_vs = 0.5,
#     )

fu = Stage(
    p_0=10.50,
    t_0=439.2,
    G_0=639,
    d_hub = 0.5,
    n = 50,
    reaction=0.235,
    alpha_0 = 90,
    alpha_1eef=14,
    H_0=38.33 - np.power(57, 2)/2000,
    c_0 = 57,
    Delta_pr=0.004,
    kappa_vs=1
)

pprint(fu.get_results())
# ff = Stage(
#     p_0=8.73,
#     t_0=570,
#     G_0=38.88,
#     d_hub=0.5,
#     n=95,
#     reaction=0.106,
#     alpha_0=90,
#     alpha_1eef=11,
#     H_0=80.39,
#     c_0=0,
#     Delta_pr=0.003,
#     extra_param={'u/c_f':0.491}
# )
# pprint(ff.get_results())

# p_vec = np.arange(0,25,1)

# for p in p_vec:
#     loc_p = ThPoint(p=p,t=550)
#     print(loc_p.k())


# def cal_time_():
#     range_ = range(500, 600, 1)
#     delta_time_vec = []
#     for H_0 in range(500, 600, 1):
#         loc_b_time = time.time()
#         f = Stage(
#             p_0=25,#1.103,
#             t_0=410.7,
#             G_0=185.7,
#             d_hub=2.275-0.133,
#             n = 50,
#             reaction=0.337,
#             alpha_0 = 90,
#             alpha_1=10,
#             H_0 = H_0,#184.5,
#             c_0 = 0,
#             Delta_pr=0.013,
#             kappa_vs = 0.5,
#             isStageOptimal= True
#         )
#         loc_e_time = time.time()
#         delta_time_vec.append(loc_e_time - loc_b_time)
#         print('time:\t',loc_e_time - loc_b_time, '\tM = ', f.M_i)
#     print('ave_time:',np.mean(delta_time_vec))
    
# cal_time_()

# #pprint(f.get_results())