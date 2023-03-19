from API.Stage import Stage
from API.TCPv2 import ThPoint

from pprint import pprint
import numpy as np

print("a")
import time

f = Stage(
    p_0=9,#1.103,
    t_0=550,
    G_0=185.7,
    d_hub=2.275-0.133,
    n = 50,
    reaction=0.337,
    alpha_0 = 90,
    alpha_1eef=10,
    H_0 = 400,#184.5,
    c_0 = 0,
    Delta_pr=0.013,
    kappa_vs = 0.5,
    isStageOptimal= True
    )

p_vec = np.arange(0,25,1)

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