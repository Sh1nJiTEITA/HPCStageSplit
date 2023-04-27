

from . import snj_solvers
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import pint
from pint import get_application_registry, Quantity
un = get_application_registry()
Q = Quantity

from datetime import datetime

import os

import pandas as pd
import numpy as np

# Термодинамические свойства воды и водяного пара
from . import TCPv2
from itertools import product

from threading import Thread

# def generate_tex_hpc(
    
# ):


GLOBAL_DF2 = pd.DataFrame(columns=[
    'n',
       
    'isOK_1',
    'is_d_k_OK_1',
    'rho_k_1',
    'd_1_1',
    'd_next_1',
    'alpha_1eef_1',
    'Z_1',
    'Z_new_1',
    'Z_ratio_1',
    'q_t_1',
    'Delta_H_1',
    'H_ave_1',
    'l_11_1',
    'd_k_1',
    'theta_1',
        
    'isOK_2',
    'is_d_k_OK_2',
    'rho_k_2',
    'd_1_2',
    'd_next_2',
    'alpha_1eef_2',
    'Z_2',
    'Z_new_2',
    'Z_ratio_2',
    'q_t_2',
    'Delta_H_2',
    'H_ave_2',
    'l_11_2',
    'd_k_2',
    'theta_2',
        
    'isOK_3',
    'is_d_k_OK_3',
    'rho_k_3',
    'd_1_3',
    'd_next_3',
    'alpha_1eef_3',
    'Z_3',
    'Z_new_3',
    'Z_ratio_3',
    'q_t_3',
    'Delta_H_3',
    'H_ave_3',
    'l_11_3',
    'd_k_3',
    'theta_3'
])

GLOBAL_DF = pd.DataFrame(columns=[
        'n', 
        'rho_k', 
                                   
        'd_1_1', 
        'd_next_1', 
        'alpha_1eef_1',
        'Z_1',
        'Z_new_1',
        'Z_ratio_1',
        'q_t_1',
        'Delta_H_1',
        'H_ave_1',
        'l_11_1',
        'd_k_1',
        'theta_1',
        
        
        'd_1_2', 
        'd_next_2', 
        'alpha_1eef_2',
        'Z_2',
        'Z_new_2',
        'Z_ratio_2',
        'q_t_2',
        'Delta_H_2',
        'H_ave_2',
        'l_11_2',
        'd_k_2',
        'theta_2',
        
        'd_1_3',
        'd_next_3', 
        'alpha_1eef_3',
        'Z_3',
        'Z_new_3',
        'Z_ratio_3',
        'q_t_3',
        'Delta_H_3',
        'H_ave_3',
        'l_11_3',
        'd_k_3',
        'theta_3'
        ])

def get_GLOBAL_DF():
    return GLOBAL_DF
def set_GLOBAL_DF(data):
    GLOBAL_DF = pd.concat([data])

class HPC_SPLIT_SETUP:
    def __init__(
        self,                       
        fp_0,                       # Давление свежего пара после гидр. сопротивления      |
        fh_0,                       # Энтальпия свежего пара                               |
        n,                          # Номинальная частота вращения                         |
        G_0,                        # Расход пара через группу ступеней                    |
        p_z,                        # Давление за группой ступеней                         |
        etaHPC_oi,                  # Внутренний КПД ЦВД                                   |
        alpha_1eef = Q(12, 'deg'),  # Эффективный угол выхода потока из сопловой решетки   |
        d_1 = Q(0.9, 'm'),          # Средней диаметр первой ступени                       |
        rho_k = Q(0.04, ''),        # Степень реактивности в корне первой ступени          |
        phi = np.mean([0.93, 0.96]),# Коэффициент скорости сопловой решетки первой ступени |
        mu = np.mean([0.95,0.97]),  # Коэффициент расхода сопловой решетки первой ступени  |
        Delta = Q(0.003, "meter"),  # Перекрыша между лопаткиами первой ступени            |
        Z = 6,                      # Преsдполагаемое число ступеней                       |
        MODE = "OVERLOAD",          # Режим расчета ступеней, перегруженный "OVERLOAD"     |
                                    # или недогруженный "UNDERLOAD"                        |
                                    #                                                      |
        isSectionFirst = False,     # Если 0, то K=[1,0.95], Если 1, то K=[0,95...]        |
        FAST_MODE = False,          # Если 0, то режим будет учитывать единицы измерения   |
        ID = 0,
        
        ):
        self.fp_0 = fp_0
        self.fh_0 = fh_0
        self.n    = n
        self.G_0  = G_0
        self.p_z  = p_z
        self.alpha_1eef = alpha_1eef
        self.etaHPC_oi = etaHPC_oi
        self.d_1 = d_1
        self.rho_k = rho_k
        self.phi = phi
        self.mu = mu
        self.Delta = Delta
        self.Z = Z
        self.MODE = MODE
        self.isSectionFirst = isSectionFirst
        self.FAST_MODE = FAST_MODE
        self.ID = ID


            
            
def local_calculate_split_smart(path,thread_number,
                                setup_1:HPC_SPLIT_SETUP, setup_2:HPC_SPLIT_SETUP, setup_3:HPC_SPLIT_SETUP, 
                                d_1_vec_1_, alpha_1eef_vec_, rho_k, 
                                delta_d_2, delta_d_3, 
                                accurancy = 2, d_step = 0.05):
    
    
    main_hub_reaction = Q(rho_k, '')
    
    setup_1.rho_k = main_hub_reaction
    setup_2.rho_k = main_hub_reaction
    setup_3.rho_k = main_hub_reaction
    
    pr = product(d_1_vec_1_,alpha_1eef_vec_)
    pr_number = int(len(list(pr)))
    
    
    max_len = int(len(list(product(d_1_vec_1_,alpha_1eef_vec_)))**4.5)
    
    main_list_1 = split_parameters_array_list(size=max_len)
    main_list_2 = split_parameters_array_list(size=max_len)
    main_list_3 = split_parameters_array_list(size=max_len)
    
    
    
    for index, param_pack_1 in enumerate(product(d_1_vec_1_, alpha_1eef_vec_)):
        SETUP_1 = setup_1
        SETUP_1.d_1 = Q(param_pack_1[0], 'm')
        SETUP_1.alpha_1eef = Q(param_pack_1[1], 'deg')
        
        SPLIT_1 = HPC_SPLIT_2(setup=SETUP_1)
        
        print('thread: {}; d_1: {}/{}'.format(thread_number, index, pr_number))
        # if split_1 is ok -> calculate split_2
        if SPLIT_1.isOK():
            # get d_k for split_1 (must be rounded)
            d_k_1 = round(SPLIT_1.get_d_k(), accurancy)
            
            # start of diameter arange for split 2 (must be rounded)
            d_start_2 = round(SPLIT_1.get_next_d(), accurancy)

            # arange of diameter values for split 2
            d_2_vec_2_ = np.arange(
                d_start_2,              # start from next_d (from split_1)
                d_start_2 + d_step + d_step * delta_d_2,  # end with previus value plus delta from input parameters
                d_step                  # which step to take between values of diameter
            )
            # print(d_2_vec_2_)
            
            for param_pack_2 in product(d_2_vec_2_, alpha_1eef_vec_):
                SETUP_2 = setup_2
                SETUP_2.d_1 = Q(param_pack_2[0], 'm')
                SETUP_2.alpha_1eef = Q(param_pack_2[1], 'deg')
        
                SPLIT_2 = HPC_SPLIT_2(setup=SETUP_2)
                # if split_2 is ok
                if SPLIT_2.isOK():
                    # get d_k for split_2 (must be rounded)
                    d_k_2 = round(SPLIT_2.get_d_k(), accurancy)
                    
                    # if diameter is increasing -> calculate split_3
                    if d_k_1 <= d_k_2:
                        
                        # diameter start of arange for split_3 (must be rounded)
                        d_start_3 = round(SPLIT_2.get_next_d(), accurancy)
                        
                        # arange of diameter values for split 3
                        d_2_vec_3_ = np.arange(
                            d_start_3,              # start from next_d (from split_2)
                            d_start_3 + d_step + d_step * delta_d_3,  # end with previus value plus delta from input parameters
                            d_step                  # which step to take between values of diameter
                        )
                        
                        for param_pack_3 in product(d_2_vec_3_, alpha_1eef_vec_):
                            SETUP_3 = setup_3
                            SETUP_3.d_1 = Q(param_pack_3[0], 'm')
                            SETUP_3.alpha_1eef = Q(param_pack_3[1], 'deg')
        
                            SPLIT_3 = HPC_SPLIT_2(setup=SETUP_3)

                            # if split_3 is ok:
                            if SPLIT_3.isOK():
                                # d_k in split_3
                                d_k_3 = round(SPLIT_3.get_d_k(), accurancy)
                                
                                # if diameter is increasing -> NICE
                                if d_k_2 <= d_k_3:
                                    main_list_1.set_row_by_split(split=SPLIT_1, isOK=1, isd_kOK=1)
                                    main_list_2.set_row_by_split(split=SPLIT_2, isOK=1, isd_kOK=1)
                                    main_list_3.set_row_by_split(split=SPLIT_3, isOK=1, isd_kOK=1)
                                
                                # if diameter is not incresing from split_2 to split_3
                                else:
                                    # continue
                                    main_list_1.set_row_by_split(split=SPLIT_1, isOK=1, isd_kOK=1)
                                    main_list_2.set_row_by_split(split=SPLIT_2, isOK=1, isd_kOK=1)
                                    main_list_3.set_row_by_split(split=SPLIT_3, isOK=1, isd_kOK=0)
                            # if split_3 is not ok:
                            else:
                                # continue
                                main_list_1.set_row_by_split(split=SPLIT_1, isOK=1, isd_kOK=1)
                                main_list_2.set_row_by_split(split=SPLIT_2, isOK=1, isd_kOK=1)
                                main_list_3.set_row_by_split(split=SPLIT_3, isOK=0, isd_kOK=-1)
                        
                    # if diameter is not increasing from split_1 to split_2:
                    else:
                        # continue 
                        main_list_1.set_row_by_split(split=SPLIT_1, isOK=1, isd_kOK=1)
                        main_list_2.set_row_by_split(split=SPLIT_2, isOK=1, isd_kOK=0)
                        main_list_3.set_row_by_split(split=SPLIT_1, isOK=-1, isd_kOK=-1)
                    
                # if split_2 is not ok:
                else:
                    # continue
                    main_list_1.set_row_by_split(split=SPLIT_1, isOK=1, isd_kOK=1)
                    main_list_2.set_row_by_split(split=SPLIT_2, isOK=0, isd_kOK=-1)
                    main_list_3.set_row_by_split(split=SPLIT_1, isOK=-1, isd_kOK=-1)
            
        # if split_1 is not ok            
        else:
            continue
            # main_list_1.set_row_by_split(split=SPLIT_1, isOK=0, isd_kOK=1)
            # main_list_2.set_row_by_split(split=SPLIT_1, isOK=-1, isd_kOK=-1)
            # main_list_3.set_row_by_split(split=SPLIT_1, isOK=-1, isd_kOK=-1)
        
    out_df = pd.DataFrame({
        'n':            setup_1.n.m, 
         
        'isOK_1':       main_list_1.isOK,
        'is_d_k_OK_1':  main_list_1.is_d_k_OK,
        'rho_k_1':      main_list_1.rho_k,
        'd_1_1':        main_list_1.d_1,
        'd_next_1':     main_list_1.d_next,
        'alpha_1eef_1': main_list_1.alpha_1eef,
        'Z_1':          main_list_1.Z,
        'Z_new_1':      main_list_1.Z_new,
        'Z_ratio_1':    main_list_1.Z_new / main_list_1.Z,
        'q_t_1':        main_list_1.q_t,
        'Delta_H_1':    main_list_1.Delta_H,
        'H_ave_1':      main_list_1.H_ave,
        'l_11_1':       main_list_1.l_11,
        'd_k_1':        main_list_1.d_k,
        'theta_1':      main_list_1.theta,
        
        'isOK_2':       main_list_2.isOK,
        'is_d_k_OK_2':  main_list_2.is_d_k_OK,
        'rho_k_2':      main_list_2.rho_k,
        'd_1_2':        main_list_2.d_1,
        'd_next_2':     main_list_2.d_next,
        'alpha_1eef_2': main_list_2.alpha_1eef,
        'Z_2':          main_list_2.Z,
        'Z_new_2':      main_list_2.Z_new,
        'Z_ratio_2':    main_list_2.Z_new / main_list_2.Z,
        'q_t_2':        main_list_2.q_t,
        'Delta_H_2':    main_list_2.Delta_H,
        'H_ave_2':      main_list_2.H_ave,
        'l_11_2':       main_list_2.l_11,
        'd_k_2':        main_list_2.d_k,
        'theta_2':      main_list_2.theta,
        
        'isOK_3':       main_list_3.isOK,
        'is_d_k_OK_3':  main_list_3.is_d_k_OK,
        'rho_k_3':      main_list_3.rho_k,
        'd_1_3':        main_list_3.d_1,
        'd_next_3':     main_list_3.d_next,
        'alpha_1eef_3': main_list_3.alpha_1eef,
        'Z_3':          main_list_3.Z,
        'Z_new_3':      main_list_3.Z_new,
        'Z_ratio_3':    main_list_3.Z_new / main_list_3.Z,
        'q_t_3':        main_list_3.q_t,
        'Delta_H_3':    main_list_3.Delta_H,
        'H_ave_3':      main_list_3.H_ave,
        'l_11_3':       main_list_3.l_11,
        'd_k_3':        main_list_3.d_k,
        'theta_3':      main_list_3.theta, 
    })
    
    out_df = out_df[out_df['rho_k_1'] != 0]
    
    out_df.to_csv(path + 'split_thread_{}.csv'.format(thread_number))


def local_calculate_split(path,thread_number,setup_1, setup_2, setup_3, d_1_vec_1_,d_1_range_3, alpha_1eef_vec_, rho_k_vec_):
    
        
                
    max_len = len(list(product(d_1_vec_1_,alpha_1eef_vec_,rho_k_vec_))) * len(alpha_1eef_vec_) * len(list(product(alpha_1eef_vec_, 
                    np.arange(0.5, d_1_range_3[1].m,d_1_range_3[2].m))))
    
    main_list_1 = split_parameters_array_list(size=max_len)
    main_list_2 = split_parameters_array_list(size=max_len)
    main_list_3 = split_parameters_array_list(size=max_len)
    
    for index, thing in enumerate(product(d_1_vec_1_,alpha_1eef_vec_,rho_k_vec_)):     #4
        setup_1_            = setup_1
        setup_1_.d_1        = Q(thing[0], 'm')
        setup_1_.alpha_1eef = Q(thing[1], 'deg')
        setup_1_.rho_k      = Q(thing[2], '')
        
        SPLIT_1 = HPC_SPLIT_2(setup=setup_1_)
        
        if (SPLIT_1.isOK()):
            for thing_ in alpha_1eef_vec_:
                setup_2_            = setup_2
                setup_2_.d_1        = Q(round(SPLIT_1.get_next_d(),3), 'm')
                setup_2_.alpha_1eef = Q(thing_, 'deg')
                setup_2_.rho_k      = Q(thing[2], '')
            
                SPLIT_2 = HPC_SPLIT_2(setup=setup_2_)
                
                if (SPLIT_2.isOK()): 
                    if (SPLIT_2.get_d_k() >= SPLIT_1.get_d_k()):
                        for index__, thing__ in enumerate(product(
                            alpha_1eef_vec_, 
                            np.arange(round(SPLIT_2.get_next_d(), 2),d_1_range_3[1].m,d_1_range_3[2].m)
                            )):
                            setup_3_            = setup_3
                            setup_3_.d_1        = Q(thing__[1], 'm') # round(SPLIT_2.get_next_d(),3) + 0
                            setup_3_.alpha_1eef = Q(thing__[0], 'deg')
                            setup_3_.rho_k      = Q(thing[2], '')
                
                            SPLIT_3 = HPC_SPLIT_2(setup=setup_3_)
                
                            if (SPLIT_3.isOK()):
                                if (SPLIT_3.get_d_k() >= SPLIT_2.get_d_k()):
                                    print('DONE')
                                    main_list_1.set_row(isOK=True, is_d_k_OK=True,
                                                        rho_k=setup_1_.rho_k.m, d_1=SPLIT_1.get_d_1(), d_next=SPLIT_1.get_next_d(), 
                                                        alpha_1eef = setup_1_.alpha_1eef.m, Z=SPLIT_1.get_Z(), Z_new = SPLIT_1.get_Z_new(), 
                                                        q_t=SPLIT_1.get_q_t(), Delta_H=SPLIT_1.get_Delta_H(), H_ave=SPLIT_1.get_H_0ave(), 
                                                        l_11=SPLIT_1.get_l_11(), d_k=SPLIT_1.get_d_k(), theta=SPLIT_1.get_last_theta())
                                    
                                    main_list_2.set_row(isOK=True, is_d_k_OK=True,
                                                        rho_k=setup_2_.rho_k.m, d_1=SPLIT_2.get_d_1(), d_next=SPLIT_2.get_next_d(), 
                                                        alpha_1eef = setup_2_.alpha_1eef.m, Z=SPLIT_2.get_Z(), Z_new = SPLIT_2.get_Z_new(), 
                                                        q_t=SPLIT_2.get_q_t(), Delta_H=SPLIT_2.get_Delta_H(), H_ave=SPLIT_2.get_H_0ave(), 
                                                        l_11=SPLIT_2.get_l_11(), d_k=SPLIT_2.get_d_k(), theta=SPLIT_2.get_last_theta())
                                    
                                    main_list_3.set_row(isOK=True, is_d_k_OK=True,
                                                        rho_k=setup_3_.rho_k.m, d_1=SPLIT_3.get_d_1(), d_next=SPLIT_3.get_next_d(), 
                                                        alpha_1eef = setup_3_.alpha_1eef.m, Z=SPLIT_3.get_Z(), Z_new = SPLIT_3.get_Z_new(), 
                                                        q_t=SPLIT_3.get_q_t(), Delta_H=SPLIT_3.get_Delta_H(), H_ave=SPLIT_3.get_H_0ave(), 
                                                        l_11=SPLIT_3.get_l_11(), d_k=SPLIT_3.get_d_k(), theta=SPLIT_3.get_last_theta())
                                else:
                                    main_list_1.set_row(isOK=True, is_d_k_OK=True,
                                                        rho_k=setup_1_.rho_k.m, d_1=SPLIT_1.get_d_1(), d_next=SPLIT_1.get_next_d(), 
                                                        alpha_1eef = setup_1_.alpha_1eef.m, Z=SPLIT_1.get_Z(), Z_new = SPLIT_1.get_Z_new(), 
                                                        q_t=SPLIT_1.get_q_t(), Delta_H=SPLIT_1.get_Delta_H(), H_ave=SPLIT_1.get_H_0ave(), 
                                                        l_11=SPLIT_1.get_l_11(), d_k=SPLIT_1.get_d_k(), theta=SPLIT_1.get_last_theta())
                                    
                                    main_list_2.set_row(isOK=True, is_d_k_OK=True,
                                                        rho_k=setup_2_.rho_k.m, d_1=SPLIT_2.get_d_1(), d_next=SPLIT_2.get_next_d(), 
                                                        alpha_1eef = setup_2_.alpha_1eef.m, Z=SPLIT_2.get_Z(), Z_new = SPLIT_2.get_Z_new(), 
                                                        q_t=SPLIT_2.get_q_t(), Delta_H=SPLIT_2.get_Delta_H(), H_ave=SPLIT_2.get_H_0ave(), 
                                                        l_11=SPLIT_2.get_l_11(), d_k=SPLIT_2.get_d_k(), theta=SPLIT_2.get_last_theta())
                                    
                                    main_list_3.set_row(isOK=True, is_d_k_OK=False,
                                                        rho_k=setup_3_.rho_k.m, d_1=SPLIT_3.get_d_1(), d_next=SPLIT_3.get_next_d(), 
                                                        alpha_1eef = setup_3_.alpha_1eef.m, Z=SPLIT_3.get_Z(), Z_new = SPLIT_3.get_Z_new(), 
                                                        q_t=SPLIT_3.get_q_t(), Delta_H=SPLIT_3.get_Delta_H(), H_ave=SPLIT_3.get_H_0ave(), 
                                                        l_11=SPLIT_3.get_l_11(), d_k=SPLIT_3.get_d_k(), theta=SPLIT_3.get_last_theta())
                                    continue   
                            else:
                                main_list_1.set_row(isOK=True, is_d_k_OK=True,
                                                    rho_k=setup_1_.rho_k.m, d_1=SPLIT_1.get_d_1(), d_next=SPLIT_1.get_next_d(), 
                                                    alpha_1eef = setup_1_.alpha_1eef.m, Z=SPLIT_1.get_Z(), Z_new = SPLIT_1.get_Z_new(), 
                                                    q_t=SPLIT_1.get_q_t(), Delta_H=SPLIT_1.get_Delta_H(), H_ave=SPLIT_1.get_H_0ave(), 
                                                    l_11=SPLIT_1.get_l_11(), d_k=SPLIT_1.get_d_k(), theta=SPLIT_1.get_last_theta())
                                    
                                main_list_2.set_row(isOK=True, is_d_k_OK=False,
                                                    rho_k=setup_2_.rho_k.m, d_1=SPLIT_2.get_d_1(), d_next=SPLIT_2.get_next_d(), 
                                                    alpha_1eef = setup_2_.alpha_1eef.m, Z=SPLIT_2.get_Z(), Z_new = SPLIT_2.get_Z_new(), 
                                                    q_t=SPLIT_2.get_q_t(), Delta_H=SPLIT_2.get_Delta_H(), H_ave=SPLIT_2.get_H_0ave(), 
                                                    l_11=SPLIT_2.get_l_11(), d_k=SPLIT_2.get_d_k(), theta=SPLIT_2.get_last_theta())
                                    
                                main_list_3.set_row(isOK=False, is_d_k_OK=-1,
                                                    d_1 = setup_3_.d_1.m,
                                                    rho_k=setup_3_.rho_k.m, alpha_1eef = setup_3_.alpha_1eef.m)
                                continue
                    else:
                        main_list_1.set_row(isOK=True, is_d_k_OK=True,
                                            rho_k=setup_1_.rho_k.m, d_1=SPLIT_1.get_d_1(), d_next=SPLIT_1.get_next_d(), 
                                            alpha_1eef = setup_1_.alpha_1eef.m, Z=SPLIT_1.get_Z(), Z_new = SPLIT_1.get_Z_new(), 
                                            q_t=SPLIT_1.get_q_t(), Delta_H=SPLIT_1.get_Delta_H(), H_ave=SPLIT_1.get_H_0ave(), 
                                            l_11=SPLIT_1.get_l_11(), d_k=SPLIT_1.get_d_k(), theta=SPLIT_1.get_last_theta())
                                    
                        main_list_2.set_row(isOK=True, is_d_k_OK=False,
                                            rho_k=setup_2_.rho_k.m, d_1=SPLIT_2.get_d_1(), d_next=SPLIT_2.get_next_d(), 
                                            alpha_1eef = setup_1_.alpha_1eef.m, Z=SPLIT_2.get_Z(), Z_new = SPLIT_2.get_Z_new(), 
                                            q_t=SPLIT_2.get_q_t(), Delta_H=SPLIT_2.get_Delta_H(), H_ave=SPLIT_2.get_H_0ave(), 
                                            l_11=SPLIT_2.get_l_11(), d_k=SPLIT_2.get_d_k(), theta=SPLIT_2.get_last_theta())
                                    
                        main_list_3.set_row(isOK=-1, is_d_k_OK=-1,
                                            rho_k=-1,#setup_3_.rho_k.m,
                                            alpha_1eef =-1)# setup_3_.alpha_1eef.m)
                        continue
                else:
                    main_list_1.set_row(isOK=True, is_d_k_OK=True,
                                        rho_k=setup_1_.rho_k.m, d_1=SPLIT_1.get_d_1(), d_next=SPLIT_1.get_next_d(), 
                                        alpha_1eef = setup_1_.alpha_1eef.m, Z=SPLIT_1.get_Z(), Z_new = SPLIT_1.get_Z_new(), 
                                        q_t=SPLIT_1.get_q_t(), Delta_H=SPLIT_1.get_Delta_H(), H_ave=SPLIT_1.get_H_0ave(), 
                                        l_11=SPLIT_1.get_l_11(), d_k=SPLIT_1.get_d_k(), theta=SPLIT_1.get_last_theta())
                                    
                    main_list_2.set_row(isOK=False, is_d_k_OK=-1,
                                        d_1 = setup_2_.d_1.m,
                                        rho_k=setup_2_.rho_k.m, 
                                        alpha_1eef = setup_2_.alpha_1eef.m)
                                    
                    main_list_3.set_row(isOK=-1, is_d_k_OK=-1,
                                        rho_k=-1,#setup_3_.rho_k.m,
                                        alpha_1eef = -1)#setup_3_.alpha_1eef.m)
                    continue    
        else:
            main_list_1.set_row(isOK=False, is_d_k_OK=True,
                                d_1 = setup_1_.d_1.m,
                                rho_k=setup_1_.rho_k.m, 
                                alpha_1eef = setup_1_.alpha_1eef.m)
                                    
            main_list_2.set_row(isOK=-1, is_d_k_OK=-1,
                                rho_k=-1,#setup_2_.rho_k.m, 
                                alpha_1eef = -1)
                                    
            main_list_3.set_row(isOK=-1, is_d_k_OK=-1,
                                rho_k=-1,#setup_3_.rho_k.m,
                                alpha_1eef = -1)# setup_3_.alpha_1eef.m)
            continue
    
    out_df = pd.DataFrame({
        'n':            setup_1_.n.m, 
         
        'isOK_1':       main_list_1.isOK,
        'is_d_k_OK_1':  main_list_1.is_d_k_OK,
        'rho_k_1':      main_list_1.rho_k,
        'd_1_1':        main_list_1.d_1,
        'd_next_1':     main_list_1.d_next,
        'alpha_1eef_1': main_list_1.alpha_1eef,
        'Z_1':          main_list_1.Z,
        'Z_new_1':      main_list_1.Z_new,
        'Z_ratio_1':    main_list_1.Z_new / main_list_1.Z,
        'q_t_1':        main_list_1.q_t,
        'Delta_H_1':    main_list_1.Delta_H,
        'H_ave_1':      main_list_1.H_ave,
        'l_11_1':       main_list_1.l_11,
        'd_k_1':        main_list_1.d_k,
        'theta_1':      main_list_1.theta,
        
        'isOK_2':       main_list_2.isOK,
        'is_d_k_OK_2':  main_list_2.is_d_k_OK,
        'rho_k_2':      main_list_2.rho_k,
        'd_1_2':        main_list_2.d_1,
        'd_next_2':     main_list_2.d_next,
        'alpha_1eef_2': main_list_2.alpha_1eef,
        'Z_2':          main_list_2.Z,
        'Z_new_2':      main_list_2.Z_new,
        'Z_ratio_2':    main_list_2.Z_new / main_list_2.Z,
        'q_t_2':        main_list_2.q_t,
        'Delta_H_2':    main_list_2.Delta_H,
        'H_ave_2':      main_list_2.H_ave,
        'l_11_2':       main_list_2.l_11,
        'd_k_2':        main_list_2.d_k,
        'theta_2':      main_list_2.theta,
        
        'isOK_3':       main_list_3.isOK,
        'is_d_k_OK_3':  main_list_3.is_d_k_OK,
        'rho_k_3':      main_list_3.rho_k,
        'd_1_3':        main_list_3.d_1,
        'd_next_3':     main_list_3.d_next,
        'alpha_1eef_3': main_list_3.alpha_1eef,
        'Z_3':          main_list_3.Z,
        'Z_new_3':      main_list_3.Z_new,
        'Z_ratio_3':    main_list_3.Z_new / main_list_3.Z,
        'q_t_3':        main_list_3.q_t,
        'Delta_H_3':    main_list_3.Delta_H,
        'H_ave_3':      main_list_3.H_ave,
        'l_11_3':       main_list_3.l_11,
        'd_k_3':        main_list_3.d_k,
        'theta_3':      main_list_3.theta, 
    })
    
    out_df = out_df[out_df['rho_k_1'] != 0]
    
    out_df.to_csv(path + 'split_thread_{}.csv'.format(thread_number))
    #print(out_df)
    
    
    #GLOBAL_DF2 = pd.concat([GLOBAL_DF2, out_df], ignore_index=True)
    
    # print(GLOBAL_DF2)
    
    #return out_df
    #GLOBAL_DF = pd.concat([local_df])
    
    #print(GLOBAL_DF)
    
    #set_GLOBAL_DF(local_df)
    
    #print(GLOBAL_DF)
    
    #return GLOBAL_DF

def HPC_SPLOT_SAME_PARAM_RANGES_FAST_SMART(
    setup_1:HPC_SPLIT_SETUP,
    setup_2:HPC_SPLIT_SETUP,
    setup_3:HPC_SPLIT_SETUP,
    thread_number,
    path = '',
    rho_k = 0.07,
    alpha_1eef_range = [9, 30, 1],
    d_start_1 = 0.5,
    delta_d_1 = 3,
    delta_d_2 = 3,
    delta_d_3 = 3,
    step = 0.05,
    accurancy = 2,
):
    
    alpha_1eef_vec = np.arange(*alpha_1eef_range)
      
    d_1_vec_1_ = np.arange(d_start_1, d_start_1 + step * delta_d_1 + step, step)
    
    d_splits = np.array_split(d_1_vec_1_, thread_number)
    
    #local_calculate_split(setup_1, setup_2, setup_3, d_splits[0], d_1_range_3, alpha_1eef_vec_, rho_k_vec_)
    
    threads = []
    
    # path,thread_number,
    # setup_1:HPC_SPLIT_SETUP, setup_2:HPC_SPLIT_SETUP, setup_3:HPC_SPLIT_SETUP, 
    # d_start_1, alpha_1eef_vec_, rho_k, 
    # delta_d_1, delta_d_2, delta_d_3, 
    # accurancy = 2, d_step = 0.5
    
    # local_calculate_split_smart(path, 0, setup_1, setup_2, setup_3,
    #                        d_splits[0], alpha_1eef_vec, rho_k,
    #                        delta_d_2, delta_d_3,
    #                        accurancy, step)
    
    for i in range(0, thread_number):
        threads.append(
            Thread(target = local_calculate_split_smart,
                   args = (path, i,
                           setup_1, setup_2, setup_3,
                           d_splits[i], alpha_1eef_vec, rho_k,
                           delta_d_2, delta_d_3,
                           accurancy, step)
            )
        )
    
    # threads start
    for i in range(0, thread_number):
        threads[i].start()
    
    # threads join
    for i in range(0, thread_number):
        threads[i].join()
    
    global GLOBAL_DF2
    
    dfs = GLOBAL_DF2
    dfs_list = []
    # combine dataframes together
    for i in range(0, thread_number):
        dfs_list.append(
            pd.read_csv(path + 'split_thread_{}.csv'.format(i))
        )
        
        
    dfs = pd.concat([dfs, *dfs_list], ignore_index=True)
    
    now = datetime.now()
    date_time = str(now.strftime("%d_%m_%H_%M_%S"))
    
    # print("path: ", )
    dfs.to_csv(path + 'split_thread_{}.csv'.format(date_time))    
    
    
    for i in range(0, thread_number):
        os.remove(path + 'split_thread_{}.csv'.format(i))
    
    return True

def HPC_SPLIT_SAME_PARAM_RANGES_FAST(
    setup_1:HPC_SPLIT_SETUP,
    setup_2:HPC_SPLIT_SETUP,
    setup_3:HPC_SPLIT_SETUP,
    path = '',
    thread_number = 3,
    d_1_range_1 = Q([0.6, 1.4, 0.1], 'm'), 
    d_1_range_3 = Q([0.6, 1.4, 0.1], 'm'),
    alpha_1eef_range = Q([9, 16, 1], 'deg'),
    rho_k_range = Q([0.03, 0.07, 0.01], "" )
):
    if path == '':return False
    
    if (isinstance(d_1_range_1,pint.Quantity)):         d_1_vec_1_ = np.arange(*(d_1_range_1.m))
    else:                                               d_1_vec_1_ = np.arange(*(d_1_range_1))

    
    
    if (isinstance(alpha_1eef_range, pint.Quantity)):   alpha_1eef_vec_ = np.arange(*(alpha_1eef_range.m))
    else:                                               alpha_1eef_vec_ = np.arange(*(alpha_1eef_range))

    if (isinstance(rho_k_range, pint.Quantity)):        rho_k_vec_ = np.arange(*(rho_k_range.m))
    else:                                               rho_k_vec_ = np.arange(*(rho_k_range))

    d_splits = np.array_split(d_1_vec_1_, thread_number)
    
    #local_calculate_split(setup_1, setup_2, setup_3, d_splits[0], d_1_range_3, alpha_1eef_vec_, rho_k_vec_)
    
    threads = []
    
    for i in range(0, thread_number):
        threads.append(
            Thread(target = local_calculate_split,
                   args = (
                       path,
                       i,
                       setup_1, 
                       setup_2, 
                       setup_3, 
                       d_splits[i], 
                       d_1_range_3, 
                       alpha_1eef_vec_, 
                       rho_k_vec_
                   )
            )
        )
    
    # threads start
    for i in range(0, thread_number):
        threads[i].start()
    
    # threads join
    for i in range(0, thread_number):
        threads[i].join()
    
    global GLOBAL_DF2
    
    dfs = GLOBAL_DF2
    dfs_list = []
    # combine dataframes together
    for i in range(0, thread_number):
        dfs_list.append(
            pd.read_csv(path + 'split_thread_{}.csv'.format(i))
        )
        
        
    dfs = pd.concat([dfs, *dfs_list], ignore_index=True)
    
    now = datetime.now()
    date_time = str(now.strftime("%d_%m_%H_%M_%S"))
    
    # print("path: ", )
    dfs.to_csv(path + 'split_thread_{}.csv'.format(date_time))    
    
    
    for i in range(0, thread_number):
        os.remove(path + 'split_thread_{}.csv'.format(i))
    
    return True
    #GLOBAL_DF = pd.concat([df1, df2, df3], ignore_index=True)
    
def HPC_SPLIT_SAME_PARAM_RANGES(
    setup_1:HPC_SPLIT_SETUP,
    setup_2:HPC_SPLIT_SETUP,
    setup_3:HPC_SPLIT_SETUP,
    d_1_range_1 = Q([0.6, 1.4, 0.1], 'm'), 
    d_1_range_3 = Q([0.6, 1.4, 0.1], 'm'),
    alpha_1eef_range = Q([9, 16, 1], 'deg'),
    rho_k_range = Q([0.03, 0.07, 0.01], "" )
    ):
    
    if (isinstance(d_1_range_1,pint.Quantity)):           d_1_vec_1 = np.arange(*(d_1_range_1.m))
    else:                                               d_1_vec_1 = np.arange(*(d_1_range_1))

    # if (isinstance(d_1_range_3,pint.Quantity)):           d_1_vec_3 = np.arange(*(d_1_range_3.m))
    # else:                                               d_1_vec_3 = np.arange(*(d_1_range_3))

    if (isinstance(alpha_1eef_range, pint.Quantity)):   alpha_1eef_vec_ = np.arange(*(alpha_1eef_range.m))
    else:                                               alpha_1eef_vec_ = np.arange(*(alpha_1eef_range))

    if (isinstance(rho_k_range, pint.Quantity)):        rho_k_vec_ = np.arange(*(rho_k_range.m))
    else:                                               rho_k_vec_ = np.arange(*(rho_k_range))

    #vec_len = len(d_1_vec_1) * len(alpha_1eef_vec_) * len(rho_k_vec_)
    
    #Z_vec           = np.zeros((vec_len))
    #Z_new_vec       = np.zeros((vec_len))
    #l_11_vec        = np.zeros((vec_len))
    #d_1_vec         = np.zeros((vec_len))
    #alpha_1eef_vec  = np.zeros((vec_len))
    #rho_k_vec       = np.zeros((vec_len))
    #next_d_vec      = np.zeros((vec_len))
    #q_t_vec         = np.zeros((vec_len))
    #Delta_H_vec     = np.zeros((vec_len))
    #____isOK_vec        = np.zeros((vec_len))
    
    out_df = pd.DataFrame(columns=[
        'n', 
        'rho_k', 
                                   
        'd_1_1', 
        'd_next_1', 
        'alpha_1eef_1',
        'Z_1',
        'Z_new_1',
        'Z_ratio_1',
        'q_t_1',
        'Delta_H_1',
        'H_ave_1',
        'l_11_1',
        'd_k_1',
        'theta_1',
        
        
        'd_1_2', 
        'd_next_2', 
        'alpha_1eef_2',
        'Z_2',
        'Z_new_2',
        'Z_ratio_2',
        'q_t_2',
        'Delta_H_2',
        'H_ave_2',
        'l_11_2',
        'd_k_2',
        'theta_2',
        
        'd_1_3',
        'd_next_3', 
        'alpha_1eef_3',
        'Z_3',
        'Z_new_3',
        'Z_ratio_3',
        'q_t_3',
        'Delta_H_3',
        'H_ave_3',
        'l_11_3',
        'd_k_3',
        'theta_3'
        ])
    
    for index, thing in enumerate(product(d_1_vec_1,         #0
                                          alpha_1eef_vec_,  #1
                                          rho_k_vec_)):     #4
        setup_1_            = setup_1
        setup_1_.d_1        = Q(thing[0], 'm')
        setup_1_.alpha_1eef = Q(thing[1], 'deg')
        setup_1_.rho_k      = Q(thing[2], '')
        
        SPLIT_1 = HPC_SPLIT_2(setup=setup_1_)
        
        if (SPLIT_1.isOK()):
            for thing_ in alpha_1eef_vec_:
                setup_2_            = setup_2
                setup_2_.d_1        = Q(round(SPLIT_1.get_next_d(),3), 'm')
                setup_2_.alpha_1eef = Q(thing_, 'deg')
                setup_2_.rho_k      = Q(thing[2], '')
            
                SPLIT_2 = HPC_SPLIT_2(setup=setup_2_)
                
                if (SPLIT_2.isOK()): 
                    if (SPLIT_2.get_d_k() >= SPLIT_1.get_d_k()):
                        for thing__ in product(
                            alpha_1eef_vec_, 
                            np.arange(round(SPLIT_2.get_next_d(), 2),d_1_range_3[1].m,d_1_range_3[2].m)
                            ):
                    
                            setup_3_            = setup_3
                            setup_3_.d_1        = Q(thing__[1], 'm') # round(SPLIT_2.get_next_d(),3) + 0
                            setup_3_.alpha_1eef = Q(thing__[0], 'deg')
                            setup_3_.rho_k      = Q(thing[2], '')
                
                            SPLIT_3 = HPC_SPLIT_2(setup=setup_3_)
                
                            if (SPLIT_3.isOK()):
                                if (SPLIT_3.get_d_k() >= SPLIT_2.get_d_k()):
                                    out_df = out_df.append({
                                    'n':            setup_1_.n.m, 
                                    'rho_k':        setup_1_.rho_k.m, 
                                   
                                    'd_1_1':        SPLIT_1.get_d_1(),
                                    'd_next_1':     SPLIT_1.get_next_d(),
                                    'alpha_1eef_1': setup_1_.alpha_1eef.m,
                                    'Z_1':          SPLIT_1.get_Z(),
                                    'Z_new_1':      SPLIT_1.get_Z_new(),
                                    'Z_ratio_1':    SPLIT_1.get_Z_new()/SPLIT_1.get_Z(),
                                    'q_t_1':        SPLIT_1.get_q_t(),
                                    'Delta_H_1':    SPLIT_1.get_Delta_H(),
                                    'H_ave_1':      SPLIT_1.get_H_0ave(),
                                    'l_11_1':       SPLIT_1.get_l_11(),
                                    'd_k_1':        SPLIT_1.get_d_k(),
                                    'theta_1':      SPLIT_1.get_last_theta(),
                                    
                                    'd_1_2':        SPLIT_2.get_d_1(), 
                                    'd_next_2':     SPLIT_2.get_next_d(),
                                    'alpha_1eef_2': setup_2_.alpha_1eef.m,
                                    'Z_2':          SPLIT_2.get_Z(),
                                    'Z_new_2':      SPLIT_2.get_Z_new(),
                                    'Z_ratio_2':    SPLIT_2.get_Z_new()/SPLIT_2.get_Z(),
                                    'q_t_2':        SPLIT_2.get_q_t(),
                                    'Delta_H_2':    SPLIT_2.get_Delta_H(),
                                    'H_ave_2':      SPLIT_2.get_H_0ave(),
                                    'l_11_2':       SPLIT_2.get_l_11(),
                                    'd_k_2':        SPLIT_2.get_d_k(),
                                    'theta_2':      SPLIT_2.get_last_theta(),

                                    'd_1_3':        SPLIT_3.get_d_1(), 
                                    'd_next_3':     SPLIT_3.get_next_d(),
                                    'alpha_1eef_3': setup_3_.alpha_1eef.m,
                                    'Z_3':          SPLIT_3.get_Z(),
                                    'Z_new_3':      SPLIT_3.get_Z_new(),
                                    'Z_ratio_3':    SPLIT_3.get_Z_new()/SPLIT_3.get_Z(),
                                    'q_t_3':        SPLIT_3.get_q_t(),
                                    'Delta_H_3':    SPLIT_3.get_Delta_H(),
                                    'H_ave_3':      SPLIT_3.get_H_0ave(),
                                    'l_11_3':       SPLIT_3.get_l_11(),
                                    'd_k_3':        SPLIT_3.get_d_k(),
                                    'theta_3':      SPLIT_3.get_last_theta()
                                },ignore_index=True)
                                    
                    else:
                        continue
                else:
                    continue
        else:
            continue
    return out_df      

def HPC_SPLIT_RANGES(                    
    fp_0 = None,                                           # Давление свежего пара после гидр. сопротивления      |
    fh_0= None,                                           # Энтальпия свежего пара                               |
    n= None,                                              # Номинальная частота вращения                         |
    G_0= None,                                            # Расход пара через группу ступеней                    |
    p_z= None,                                            # Давление за группой ступеней                         |
    etaHPC_oi= None,                                      # Внутренний КПД ЦВД                                   |
    ID= None,                                             # ...                                                  |
    alpha_1eef_range    = Q([9, 20, 1], "deg"),     # Эффективный угол выхода потока из сопловой решетки   |
    d_1_range           = Q([0.6, 1.4, 0.1], "m"),  # Средней диаметр первой ступени                       |
    rho_k_range         = Q([0.03, 0.07, 0.01], ""),# Степень реактивности в корне первой ступени          |
    phi                 = np.mean([0.93, 0.96]),    # Коэффициент скорости сопловой решетки первой ступени |
    mu                  = np.mean([0.95,0.97]),     # Коэффициент расхода сопловой решетки первой ступени  |
    Delta               = Q(0.003, "meter"),        # Перекрыша между лопаткиами первой ступени            |
                                                    # Предполагаемое число ступеней                        |
    MODE                = "OVERLOAD",               # Режим расчета ступеней, перегруженный "OVERLOAD"     |
                                                    # или недогруженный "UNDERLOAD"                        |
    isSectionFirst      = False,                    # Режим расчета первого отсека или не первого          |
    FAST_MODE           = True,                     # Режим быстрых вычислений                             |
    
    setup:HPC_SPLIT_SETUP = None
    )-> pd.DataFrame: 
    '''
        Данная функция предназначена для расчета разбивок группы ступеней ЦВД с разными начальными значениями
        (диапазонами) диаметров первой нерегулируемой ступени, эффективного угла входа в сопловую решетку, степени
        Раекивности в корне. 
        Расчет осуществляется из предположения о том, что объем линейно изменяется от ступени к ступени.
        Возвращаемым типом данных 
        является pandas.DataFrame, который включает в себя следующие столбцы:
        
        n           [1/s]   - Частота вращения вала
        d_1         [m]     - Средний диаметр первой нерегулируемой ступени
        alpha_1eef  [deg]   - Эффективный угол входа потока в сопловую решетку первой нерегулируемой ступени
        rho_k       []      - Степень реактивности в корне (одинакова по всей проточной части)
        Z           []      - Предполагаемое целое число ступеней
        Z_new       []      - Конечное неокругленное число ступеней
        Z_ratio     []      - Невязка конечного количества ступеней по отношению к предполагаемому
        q_t         []      - Коэффициент возврата теплоты
        Delta_H     [kJ/kg] - Невязка теплоперепадов, располагаемая на целые ступени
        l_11        [m]     - Высота сопловых лопаток первой нерегулируемой ступени
        d_next      [m]     - Предполагаемый средний диаметр первой ступени следующей группы ступеней
        ID          []      - ...
        isOK        []      - Возможен ли расчет ступеней

        У функции присутствуют три режима работы.
        
        1. Какую разбивку необходим посчитать:
        MODE = "OVERLOAD"   - Расчет производится для получения перегруженных ступеней
        MODE = "UNDERLOAD"  - Расчет производится для получения недогруженных ступеней
        
        2. Как разбивку необходимо посчитать:
        FAST_MODE = True    - Расчет производится без учета единиц измерения, что увеличивает скорость в 10.14 раз
        FAST_MODE = False   - Расчет проивзодится с учетом единиц измерения 
        Целесообразно производить расчет в режиме FAST_MODE = True, так как результаты получаются одинаковые, а
        скорость расчета многократно возрастает
        
        3. isSectionFirst = True - в этом случае коэффициент К начинается с 1, то есть K = [1, 0.95 ...]
        Иначе К = [0.95, 0.95 ...]
        
    '''
    
    if (isinstance(d_1_range,pint.Quantity)):           d_1_vec_ = np.arange(*(d_1_range.m))
    else:                                               d_1_vec_ = np.arange(*(d_1_range))

    if (isinstance(alpha_1eef_range, pint.Quantity)):   alpha_1eef_vec_ = np.arange(*(alpha_1eef_range.m))
    else:                                               alpha_1eef_vec_ = np.arange(*(alpha_1eef_range))

    if (isinstance(rho_k_range, pint.Quantity)):        rho_k_vec_ = np.arange(*(rho_k_range.m))
    else:                                               rho_k_vec_ = np.arange(*(rho_k_range))

    vec_len = len(d_1_vec_) * len(alpha_1eef_vec_) * len(rho_k_vec_)
    
    if not(setup == None): n_ = setup.n.m
    else: n_ = n.m
    
    Z_vec           = np.zeros((vec_len))
    Z_new_vec       = np.zeros((vec_len))
    l_11_vec        = np.zeros((vec_len))
    d_1_vec         = np.zeros((vec_len))
    alpha_1eef_vec  = np.zeros((vec_len))
    rho_k_vec       = np.zeros((vec_len))
    next_d_vec      = np.zeros((vec_len))
    q_t_vec         = np.zeros((vec_len))
    Delta_H_vec     = np.zeros((vec_len))
    d_k_vec         = np.zeros((vec_len))
    isOK_vec        = np.zeros((vec_len))
    
    setup_ = setup
    # Производим разбивку для каждой из вариаций 
    for i,it in enumerate(product(d_1_vec_,alpha_1eef_vec_,rho_k_vec_)):
        if not(setup == None):
            setup_.d_1          = Q(it[0], 'm')
            setup_.alpha_1eef   = Q(it[1], 'deg')
            setup_.rho_k        = Q(it[2], '')
            local_split = HPC_SPLIT_2(setup=setup_)
            
        else:
            local_split = HPC_SPLIT_2(
                fp_0                    =fp_0,
                fh_0                    =fh_0,
                n                       =n,
                G_0                     =G_0,
                p_z                     =p_z,
                MODE                    =MODE,
                etaHPC_oi               =etaHPC_oi,
                d_1                     =Q(it[0],'m'),
                alpha_1eef              =Q(it[1], 'deg'),
                rho_k                   =Q(it[2], ''),
                Delta                   =Delta,
                phi                     =phi,
                mu                      =mu,
                isSectionFirst          =isSectionFirst,
                FAST_MODE               =FAST_MODE      
            )
        
        # Если разбивка была просчитана 
        if (local_split.isOK()):
            # Так как в быстром режиме не ведется учет единиц измерения, то просто получаем их
            if (FAST_MODE):
                Z_vec[i]            = local_split.get_Z()
                Z_new_vec[i]        = local_split.get_Z_new()
                l_11_vec[i]         = local_split.get_l_11()
                next_d_vec[i]       = local_split.get_next_d()
                q_t_vec[i]          = local_split.get_q_t()
                Delta_H_vec[i]      = local_split.get_Delta_H()
                d_k_vec[i]          = local_split.get_d_k()
            # В обычном режиме получаем значение (без единиц измерения)
            else:
                Z_vec[i]            = local_split.get_Z()
                Z_new_vec[i]        = local_split.get_Z_new().m
                l_11_vec[i]         = local_split.get_l_11().m
                next_d_vec[i]       = local_split.get_next_d().m
                q_t_vec[i]          = local_split.get_q_t().m
                Delta_H_vec[i]      = local_split.get_Delta_H().m
                d_k_vec[i]          = local_split.get_d_k().m
            
            isOK_vec[i] = True
        # Если разбивка не получилась
        else:
            Z_vec[i]            = local_split.get_Z()
            Z_new_vec[i]        = 0
            l_11_vec[i]         = 0
            next_d_vec[i]       = 0
            q_t_vec[i]          = 0
            Delta_H_vec[i]      = 0
            d_k_vec[i]          = 0
            
            isOK_vec[i] = False
        
        d_1_vec[i]          = it[0]
        alpha_1eef_vec[i]   = it[1]
        rho_k_vec[i]        = it[2]
    
    # выходная таблица
    data = {
        'n':np.full(vec_len, n_),
        'd_1':d_1_vec,
        'alpha_1eef':alpha_1eef_vec,
        'rho_k':rho_k_vec,
        'Z':Z_vec,
        'Z_new':Z_new_vec,
        'Z_ratio':(Z_new_vec/Z_vec),
        'q_t':q_t_vec,
        'Delta_H':Delta_H_vec,
        'l_11':l_11_vec,
        'd_next':next_d_vec,
        'd_k':d_k_vec,
        'ID':np.full(vec_len, ID),
        '__isOK':isOK_vec
    }
    return pd.DataFrame(data)

class HPC_SPLIT_2:
    def __init__(
        self,                       
        fp_0 = None,                       # Давление свежего пара после гидр. сопротивления      |
        fh_0 = None,                       # Энтальпия свежего пара                               |
        n= None,                          # Номинальная частота вращения                         |
        G_0= None,                        # Расход пара через группу ступеней                    |
        p_z= None,                        # Давление за группой ступеней                         |
        alpha_1eef= None,                 # Эффективный угол выхода потока из сопловой решетки   |
        etaHPC_oi= None,                  # Внутренний КПД ЦВД                                   |
        d_1= None,                        # Средней диаметр первой ступени                       |
        rho_k= None,                      # Степень реактивности в корне первой ступени          |
        phi = np.mean([0.93, 0.96]),# Коэффициент скорости сопловой решетки первой ступени |
        mu = np.mean([0.95,0.97]),  # Коэффициент расхода сопловой решетки первой ступени  |
        Delta = Q(0.003, "meter"),  # Перекрыша между лопаткиами первой ступени            |
        Z = 6,                      # Преsдполагаемое число ступеней                       |
        MODE = "OVERLOAD",          # Режим расчета ступеней, перегруженный "OVERLOAD"     |
                                    # или недогруженный "UNDERLOAD"                        |
                                    #                                                      |
        isSectionFirst = False,     # Если 0, то K=[1,0.95], Если 1, то K=[0,95...]        |
        FAST_MODE = False,          # Если 0, то режим будет учитывать единицы измерения   |
        setup = None
    ):
        if not(setup == None):
            # Размерные величины
            self.__fp_0 = setup.fp_0
            self.__fh_0 = setup.fh_0
            self.__n    = setup.n
            self.__G_0  = setup.G_0
            self.__p_z  = setup.p_z
            self.__d_1  = setup.d_1
            self.__alpha_1eef = setup.alpha_1eef
        
            # доп
            self.__FAST_MODE  = setup.FAST_MODE
            self.__isSectionFirst = setup.isSectionFirst
        
            # Безразмерные величины
            self.__Z    = setup.Z
            self.__etaHPC_oi = setup.etaHPC_oi
            self.__rho_k = setup.rho_k
            self.__phi  = setup.phi
            self.__mu = setup.mu
            self.__Delta = setup.Delta
            self.__MODE = setup.MODE
        else:
            # Размерные величины
            self.__fp_0 = fp_0
            self.__fh_0 = fh_0
            self.__n    = n
            self.__G_0  = G_0
            self.__p_z  = p_z
            self.__d_1  = d_1
            self.__alpha_1eef = alpha_1eef
        
            # доп
            self.__FAST_MODE  = FAST_MODE
            self.__isSectionFirst = isSectionFirst
        
            # Безразмерные величины
            self.__Z    = Z
            self.__etaHPC_oi = etaHPC_oi
            self.__rho_k = rho_k
            self.__phi  = phi
            self.__mu = mu
            self.__Delta = Delta
            self.__MODE = MODE    
        
        # Out values
        self.__l_11                 = 0
        self.__d_vec                = None
        self.__l_vec                = None
        self.__stages_number_vec    = None
        self.__theta_vec            = None
        self.__rho_vec              = None
        self.__K_vec                = None
        self.__H_vec                = None
        self.__H_0ave               = None
        self.__q_t                  = None
        self.__Z_new                = None
        self.__Delta_H              = None
        self.__H_new_vec            = None
        self.__next_d               = None
        self.__d_k                  = None
        
        # Если в процессе расчета высоты первой лопатки произошло исключение - вариант не считается
        self.__isOK                   = True
        self.__isDim                = True
                  
        
        if (self.__FAST_MODE):
            self.__calculate_faster()
        else:
            self.__calculate()

    def isOK(self):
        
        # Если единицы измерения отключены
        if not(self.__isDim):
            if (self.__isOK):
                if ((self.__Z_new/self.__Z <= 1.08) & (self.__l_11 >= 0.013)):      
                        return True
                else:   return False
            else:   return False
        # Если единицы измерения включены
        else:
            if (self.__isOK):
                if ((self.__Z_new/self.__Z <= 1.08) & (self.__l_11.m >= 0.013)):      
                    return True
                else:   return False
            else: return False    
            
        
    
    def __str__(self) -> str:
        # n	d_1	alpha_1eef	rho_k	Z	Z_new	Z_ratio	q_t	Delta_H	l_11	d_next	ID	____isOK
        
        out_str  = ''
        out_str += 'n = {} \n'.format(self.__n)
        out_str += 'd_1 = {} \n'.format(self.__d_1)
        out_str += 'alpha_1eef = {} \n'.format(self.__alpha_1eef)
        out_str += 'rho_k = {} \n'.format(self.__rho_k)
        out_str += "Z = {} \n".format(self.__Z)
        out_str += "Z_new = {} \n".format(self.__Z_new)
        out_str += 'Z_ratio = {} \n'.format(self.__Z_new/self.__Z)
        out_str += "q_t = {} \n".format(self.__q_t)
        out_str += 'Delta_H = {} \n'.format(self.__Delta_H)
        out_str += "l_11 = {} \n".format(self.__l_11)
        out_str += 'd_next = {} \n'.format(self.__next_d)
        out_str += '____isOK = {} \n'.format(self.__isOK)
        out_str += 'FAST_MODE = {}'.format(self.__FAST_MODE)
        
        return out_str
        
    def __normolize_dimensions(self):
        '''
            Функция переводит внутренние заданные величины в нужные для расчета единицы измерения
            для более быстрого расчета
        '''   
        if isinstance(self.__fp_0, pint.Quantity): self.__fp_0.ito("MPa")
        elif not isinstance(self.__fp_0, int) or not isinstance(self.__fp_0, float):
            raise Exception("Input fp_0 is not a number")
        
        if isinstance(self.__fh_0, pint.Quantity): self.__fh_0.ito("kJ/kg")
        elif not isinstance(self.__fh_0, int) or not isinstance(self.__fh_0, float):
            raise Exception("Input fh_0 is not a number") 
        
        if isinstance(self.__p_z, pint.Quantity): self.__p_z.ito("MPa")
        elif not isinstance(self.__p_z, int) or not isinstance(self.__p_z, float):
            raise Exception("Input p_z is not a number") 
        
        if isinstance(self.__d_1, pint.Quantity): self.__d_1.ito("m")
        elif not isinstance(self.__d_1, int) or not isinstance(self.__d_1, float):
            raise Exception("Input d_1 is not a number") 
        
        if isinstance(self.__alpha_1eef, pint.Quantity): self.__alpha_1eef.ito("deg")
        elif not isinstance(self.__alpha_1eef, int) or not isinstance(self.__alpha_1eef, float):
            raise Exception("Input alpha_1eef is not a number")
        
        if isinstance(self.__n, pint.Quantity): self.__n.ito("1/s")
        elif not isinstance(self.__n, int) or not isinstance(self.__n, float):
            raise Exception("Input n is not a number")    
        
    def __turn_off_input_parameters_dimension(self):
        ""
        self.__fp_0         = self.__fp_0.m
        self.__fh_0         = self.__fh_0.m
        self.__n            = self.__n.m
        self.__G_0          = self.__G_0.m
        self.__p_z          = self.__p_z.m
        self.__d_1          = self.__d_1.m
        self.__alpha_1eef   = self.__alpha_1eef.m
        self.__rho_k        = self.__rho_k.m
        self.__Delta        = self.__Delta.m
        self.__isDim        = False
    
    # TODO
    def __turn_on_input_parameters_dimension(self):
        self.__fp_0         = Q(self.__fp_0,        'MPa')
        self.__fh_0         = Q(self.__fh_0,        'kJ/kg')
        self.__n            = Q(self.__n,           '1/s')
        self.__G_0          = Q(self.__G_0,         'kg/s')
        self.__p_z          = Q(self.__p_z,         'MPa')
        self.__d_1          = Q(self.__d_1,         'm')
        self.__alpha_1eef   = Q(self.__alpha_1eef,  'deg')
        self.__rho_k        = Q(self.__rho_k,       '')
        self.__Delta        = Q(self.__Delta,       'm')
        self.__isDim        = True
    
    def __calculate_faster(self):
        try:
            if (self.__isDim): self.__turn_off_input_parameters_dimension()
            
            if (self.__Z == 0): raise Exception('z =< 0')
            if (self.__Z == 1): raise Exception('z == 1')
            # 3.1
            # ...
            
            # 3.2
            fPoint_0    = TCPv2.ThPoint(p=self.__fp_0, h=self.__fh_0)
        
            # 3.3
            theta = 20
            iterations_number = 0
            percent_difference = 0.0
            
            while (True):
    
                # 3.2. Определяем степень реактивности на среднем диаметре
                rho = self.__rho_k + 1.8/(theta+1.8)

                if (rho > 1): raise Exception("RHO > 1")
                # 3.3. Определяем оптимальное значение u/c_f
                uDIVu_cf = self.__phi*np.cos(np.radians(self.__alpha_1eef))/(2*np.sqrt(1-rho))
    
                # 3.4. Определяем располагаемый теплопеперад по параметрам торможения при оптимальном
                # u/c_f для первой ступени 
                fH_01 = 12300 * np.power((self.__d_1 * self.__n)/(uDIVu_cf * 50), 2)/1000
                #fH_01.ito('kJ/kg')
    
                # 3.5. Определяем теоретическое значение энтальпии за первой ступенью
                h_2t = self.__fh_0 - fH_01
    
                # 3.6. Определяем удальный объем пара за первой нерегулируемой ступенью при
                # изоэнтропном процессе расширения по свойствам воды и водяного пара
                Point_2t = TCPv2.ThPoint(h=h_2t, s=fPoint_0.s(0))
                v_2t = Point_2t.v(0)
    
                # 3.7. Определеяем высоту первой нерегулируемой ступени
                self.__l_11 = (self.__G_0 * v_2t * uDIVu_cf)/(np.power(np.pi*self.__d_1,2) * self.__n * np.sqrt(1-rho) * np.sin(np.radians(self.__alpha_1eef)) * self.__mu)
                
                # 3.8. Определяем окончательное значение обратной веерности и проверяем его
    
                # Проверка условия 
                # Если получившаяся величина больше, чем заданная
                if (self.__d_1/self.__l_11 > theta):
                    # То вычитаем отношение из единицы
                    percent_difference = (1 - theta/(self.__d_1/self.__l_11)) 
                # Если величина меньше, чем заданная
                else:
                    # То вычитаем из отношения единицу
                    percent_difference = (theta/(self.__d_1/self.__l_11) - 1)
    
                # Если условие выполнилось - выходим из цикла
                if (np.abs(percent_difference) < 0.01):
                    break
                # Иначе добавляем итерацию и меняем приближающее значение на найденное в процессе цикла
                else:
                    iterations_number += 1
                    theta = (self.__d_1/self.__l_11)   
            
            # 4.
            # Высота рабочей лопатки
            l_21 = self.__l_11 + self.__Delta
            #l_21.ito('m')
            
            # 5.
            # Корневой диаметр ступени
            self.__d_k = self.__d_1 - l_21
            #d_k.ito('m')
            
            # 6.1. Значение энтальпии пара при изоэнтропном расширении пара в ЦВД:
            # Термодинамическая точка zt
            Point_zt = TCPv2.ThPoint(p=self.__p_z, s=fPoint_0.s(0))
            # Энтальпия
            h_zt = Point_zt.h(0)

            # 6.2. Теоретический перепад на отсек нерегулируемых ступеней ЦВД:
            fH_0 = self.__fh_0 - h_zt

            # 6.3. Действительный теплоперепад на отсек нерегулируемых ступеней ЦВД
            H_i = fH_0 * self.__etaHPC_oi

            # 6.4. Действительное значение энтальпии за ЦВД (за последней ступенью)
            h_z = self.__fh_0 - H_i

            # 6.5 Действительный объем за ЦВД (за последней ступенью)
            # Термодинамическая точка 2z
            Point_2z = TCPv2.ThPoint(p=self.__p_z, h=h_z)
            v_2z = Point_2z.v(0)
            
            # 7 Высота рабочей лопатки последней ступени
            l_2z = (-self.__d_k + np.sqrt(self.__d_k**2 + 4*l_21*self.__d_1*v_2z/v_2t))/2
            
            # 8 Средний диаметр последней ступени группы
            d_2z = self.__d_k + l_2z
            
            # 9.1 Обратная вверность в первой и последней ступени
            #theta_1 = (l_21 + d_k)/l_21
            #theta_z = (l_2z + d_k)/l_2z
            # 9.2 Степень реактивности в первой и последней ступени
            #rho_1 = self.__rho_k + 1.8/(theta_1 + 1.8)
            #rho_z = self.__rho_k + 1.8/(theta_z + 1.8)
            # 9.3 Оптимальное значение u/c_ф
            #uDIVu_1 = self.__phi*np.cos(self.__alpha_1eef)/(2*np.sqrt(1-rho_1))
            #uDIVu_z = self.__phi*np.cos(self.__alpha_1eef)/(2*np.sqrt(1-rho_z))
            
            ## 10.1.1, 10.1.2
            # Коэффициент k для функции диаметра
            # k_d = (d_2z - self.__d_1)/(self.__Z-1)
            # b_d = (self.__d_1*self.__Z - d_2z)/(self.__Z-1)

            # Коэффициент k для функции длины лопатки
            k_l = (l_2z - l_21)/(self.__Z-1)
            b_l = (l_21*self.__Z - l_2z)/(self.__Z-1)

            # Функция для определения диаметра по номеру ступени
            # def d(z_i):
            #     return k_d*z_i + b_d

            # Функция для определения длины лопатки по номеру ступени
            def l(z_i):
                return k_l*z_i + b_l

            # Перечислим номера ступеней в данном векторе
            self.__stages_number_vec = np.arange(1, self.__Z+1, 1)

            # Длины сопловых лопаток каждой ступени
            self.__l_vec = l(self.__stages_number_vec)
            
            # Диаметры каждой ступени
            self.__d_vec = self.__l_vec + self.__d_k

            # 10.2
            # Обратная вверность для каждой ступени
            self.__theta_vec = (self.__l_vec + self.__d_k)/self.__l_vec

            # 10.3
            # Степень реактивности на среднем диаметре для каждой ступени
            self.__rho_vec = self.__rho_k + 1.8/(self.__theta_vec + 1.8)
            
            # 10.4. Для каждой ступени определяем величину u/c_f
            uDIVc_f_vec = self.__phi*np.cos(np.radians(self.__alpha_1eef))/(2*np.sqrt(1-self.__rho_vec))

            # 10.5 Теплоперепад по статическим параметрам для каждой ступени
            # Вектор коэффициентов K_i
            # k
            
            self.__K_vec = np.full(self.__Z, 0.95)
            
            if (self.__isSectionFirst): self.__K_vec[0] = 1.0
                
            self.__H_vec = 12300 * (self.__d_vec/uDIVc_f_vec)**2 * (self.__n/50)**2 * self.__K_vec/1000
            #self.__H_vec.ito('kJ/kg')
            
            # 10.6 Среднее значение теплоперепада за группу ступеней
            self.__H_0ave = np.mean(self.__H_vec)
            
            # 10.7 Коэффициент возврата теплоты
            self.__q_t = 4.8*10**(-4) * (1 - self.__etaHPC_oi)*self.__H_0ave * (self.__Z-1)/self.__Z
            
            # 10.8 Уточненное количество ступеней группы
            self.__Z_new = fH_0/self.__H_0ave * (1+self.__q_t)
            
            #self.__Z_new.ito('')

            # 11 Величина распределения теплоперепадов дробной ступени на остальные
            self.__Delta_H = (fH_0*(1+self.__q_t)/self.__Z) - self.__H_0ave
            
            # 12 Уточненные теплоперепады
            #self.__H_new_vec = self.__H_vec + self.__Delta_H
            
            # Диаметр первой ступени следующего цилиндра
            self.__next_d = self.__d_vec[len(self.__d_vec)-1] + (self.__d_vec[1]-self.__d_vec[0])
            
            # Для автоматическогов выбора ступеней
            if ((self.__MODE == "OVERLOAD")and not(abs(int(self.__Z_new) - self.__Z) < 1)):
                self.__Z = int(self.__Z_new)
                self.__calculate_faster()

            if ((self.__MODE == "UNDERLOAD")and not(int(self.__Z_new) - self.__Z > 1)):
                self.__Z = int(self.__Z_new)+1
                self.__calculate_faster()
                      
        except Exception as s:
            self.__isOK = False
            #print(s)     
        
    def __calculate(self):
        try:
            # 3.1
            # ...
            
            # 3.2
            fPoint_0    = TCPv2.ThPoint(p=self.__fp_0, h=self.__fh_0)
        
            # 3.3
            theta = 20
            iterations_number = 0
            percent_difference = 0.0
            
            while (True):
    
                # 3.2. Определяем степень реактивности на среднем диаметре
                rho = self.__rho_k + 1.8/(theta+1.8)
    
                # 3.3. Определяем оптимальное значение u/c_f
                uDIVu_cf = self.__phi*np.cos(self.__alpha_1eef)/(2*np.sqrt(1-rho))
    
                # 3.4. Определяем располагаемый теплопеперад по параметрам торможения при оптимальном
                # u/c_f для первой ступени 
                fH_01 = 12300 * np.power((self.__d_1 * self.__n)/(uDIVu_cf * 50), 2)
                fH_01.ito('kJ/kg')
    
                # 3.5. Определяем теоретическое значение энтальпии за первой ступенью
                h_2t = self.__fh_0 - fH_01
    
                # 3.6. Определяем удальный объем пара за первой нерегулируемой ступенью при
                # изоэнтропном процессе расширения по свойствам воды и водяного пара
                Point_2t = TCPv2.ThPoint(h=h_2t, s=fPoint_0.s())
                v_2t = Point_2t.v()
    
                # 3.7. Определеяем высоту первой нерегулируемой ступени
                self.__l_11 = (self.__G_0 * v_2t * uDIVu_cf)/(np.power(np.pi*self.__d_1,2) * self.__n * np.sqrt(1-rho) * np.sin(self.__alpha_1eef) * self.__mu)
                
                # 3.8. Определяем окончательное значение обратной веерности и проверяем его
    
                # Проверка условия 
                # Если получившаяся величина больше, чем заданная
                if (self.__d_1/self.__l_11 > theta):
                    # То вычитаем отношение из единицы
                    percent_difference = (1 - theta/(self.__d_1/self.__l_11)) 
                # Если величина меньше, чем заданная
                else:
                    # То вычитаем из отношения единицу
                    percent_difference = (theta/(self.__d_1/self.__l_11) - 1)
    
                # Если условие выполнилось - выходим из цикла
                if (np.abs(percent_difference) < 0.01):
                    break
                # Иначе добавляем итерацию и меняем приближающее значение на найденное в процессе цикла
                else:
                    iterations_number += 1
                    theta = (self.__d_1/self.__l_11)
            
            self.__l_11.ito("m")   
            
            # 4.
            # Высота рабочей лопатки
            l_21 = self.__l_11 + self.__Delta
            l_21.ito('m')
            
            # 5.
            # Корневой диаметр ступени
            self.__d_k = self.__d_1 - l_21
            self.__d_k.ito('m')
            
            # 6.1. Значение энтальпии пара при изоэнтропном расширении пара в ЦВД:
            # Термодинамическая точка zt
            Point_zt = TCPv2.ThPoint(p=self.__p_z, s=fPoint_0.s())
            # Энтальпия
            h_zt = Point_zt.h()

            # 6.2. Теоретический перепад на отсек нерегулируемых ступеней ЦВД:
            fH_0 = self.__fh_0 - h_zt

            # 6.3. Действительный теплоперепад на отсек нерегулируемых ступеней ЦВД
            H_i = fH_0 * self.__etaHPC_oi

            # 6.4. Действительное значение энтальпии за ЦВД (за последней ступенью)
            h_z = self.__fh_0 - H_i

            # 6.5 Действительный объем за ЦВД (за последней ступенью)
            # Термодинамическая точка 2z
            Point_2z = TCPv2.ThPoint(p=self.__p_z, h=h_z)
            v_2z = Point_2z.v()
            
            # 7 Высота рабочей лопатки последней ступени
            l_2z = (-self.__d_k + np.sqrt(self.__d_k**2 + 4*l_21*self.__d_1*v_2z/v_2t))/2
            
            # 8 Средний диаметр последней ступени группы
            d_2z = self.__d_k + l_2z
            
            # 9.1 Обратная вверность в первой и последней ступени
            theta_1 = (l_21 + self.__d_k)/l_21
            theta_z = (l_2z + self.__d_k)/l_2z
            # 9.2 Степень реактивности в первой и последней ступени
            rho_1 = self.__rho_k + 1.8/(theta_1 + 1.8)
            rho_z = self.__rho_k + 1.8/(theta_z + 1.8)
            # 9.3 Оптимальное значение u/c_ф
            uDIVu_1 = self.__phi*np.cos(self.__alpha_1eef)/(2*np.sqrt(1-rho_1))
            uDIVu_z = self.__phi*np.cos(self.__alpha_1eef)/(2*np.sqrt(1-rho_z))
            
            ## 10.1.1, 10.1.2
            # Коэффициент k для функции диаметра
            k_d = (d_2z - self.__d_1)/(self.__Z-1)
            b_d = (self.__d_1*self.__Z - d_2z)/(self.__Z-1)

            # Коэффициент k для функции длины лопатки
            k_l = (l_2z - l_21)/(self.__Z-1)
            b_l = (l_21*self.__Z - l_2z)/(self.__Z-1)

            # Функция для определения диаметра по номеру ступени
            def d(z_i):
                return k_d*z_i + b_d

            # Функция для определения длины лопатки по номеру ступени
            def l(z_i):
                return k_l*z_i + b_l

            # Перечислим номера ступеней в данном векторе
            self.__stages_number_vec = np.arange(1, self.__Z+1, 1)

            # Диаметры каждой ступени
            self.__d_vec = d(self.__stages_number_vec)

            # Длины сопловых лопаток каждой ступени
            self.__l_vec = l(self.__stages_number_vec)
            
            # 10.2
            # Обратная вверность для каждой ступени
            self.__theta_vec = (self.__l_vec + self.__d_k)/self.__l_vec

            # 10.3
            # Степень реактивности на среднем диаметре для каждой ступени
            self.__rho_vec = self.__rho_k + 1.8/(self.__theta_vec + 1.8)
            
            # 10.4. Для каждой ступени определяем величину u/c_f
            uDIVc_f_vec = self.__phi*np.cos(self.__alpha_1eef)/(2*np.sqrt(1-self.__rho_vec))

            # 10.5 Теплоперепад по статическим параметрам для каждой ступени
            # Вектор коэффициентов K_i
            # k
            
            self.__K_vec = np.full(self.__Z, 0.95)
            
            if (self.__isSectionFirst): self.__K_vec[0] = 1.0
                
            self.__H_vec = 12300 * (self.__d_vec/uDIVc_f_vec)**2 * (self.__n/50)**2 * self.__K_vec
            self.__H_vec.ito('kJ/kg')
            
            # 10.6 Среднее значение теплоперепада за группу ступеней
            self.__H_0ave = np.mean(self.__H_vec)
            
            # 10.7 Коэффициент возврата теплоты
            q_t_k = Q(4.8*10**(-4), 'kg/kJ')
            self.__q_t = q_t_k * (1 - self.__etaHPC_oi)*self.__H_0ave * (self.__Z-1)/self.__Z
            
            # 10.8 Уточненное количество ступеней группы
            self.__Z_new = fH_0/self.__H_0ave * (1+self.__q_t)
            self.__Z_new.ito('')

            # 11 Величина распределения теплоперепадов дробной ступени на остальные
            self.__Delta_H = (fH_0*(1+self.__q_t)/self.__Z) - self.__H_0ave
            
            # 12 Уточненные теплоперепады
            self.__H_new_vec = self.__H_vec + self.__Delta_H
            
            # Диаметр первой ступени следующего цилиндра
            self.__next_d = self.__d_vec[len(self.__d_vec)-1] + (self.__d_vec[1]-self.__d_vec[0])
            
            # Для автоматическогов выбора ступеней
            if ((self.__MODE == "OVERLOAD")and not(abs(int(self.__Z_new.m) - self.__Z) < 1)):
                self.__Z = Q(int(self.__Z_new), "").m
                self.__calculate()
            #elif (((self.__MODE == "OVERLOAD")and (int(self.__Z_new.m) == self.__Z + 1))):
                #pass# not(int(self.__Z_new.m)+1 == self.__Z)
            if ((self.__MODE == "UNDERLOAD")and not(int(self.__Z_new.m) - self.__Z > 1)):
                self.__Z = Q(int(self.__Z_new) + 1, "").m
                self.__calculate()
            #else:
                #pass
                
            
        except Exception as s:
            self.__isOK = False
            #print(s)
    
    # getters
    def get_last_theta(self):           return self.__theta_vec[len(self.__theta_vec) - 1]
    def get_d_k(self):                  return self.__d_k
    def get_d_1(self):                  return self.__d_1
    def get_q_t(self):                  return self.__q_t
    def get_Delta(self):                return self.__Delta_H
    def get_next_d(self):               return self.__next_d
    def get_Z(self):                    return self.__Z
    def get_l_11(self):                 return self.__l_11
    def get_d_vec(self):                return self.__d_vec
    def get_l_vec(self):                return self.__l_vec
    def get_stages_number_vec(self):    return self.__stages_number_vec
    def get_theta_vec(self):            return self.__theta_vec
    def get_rho_vec(self):              return self.__rho_vec
    def get_K_vec(self):                return self.__K_vec
    def get_H_vec(self):                return self.__H_vec
    def get_H_0ave(self):               return self.__H_0ave
    def get_q_t(self):                  return self.__q_t          
    def get_Z_new(self):                return self.__Z_new          
    def get_Delta_H(self):              return self.__Delta_H       
    def get_H_new_vec(self):            return self.__H_new_vec      
    def get_rho_k(self):                return self.__rho_k
    def get_alpha_1eef(self):           return self.__alpha_1eef
    
class split_parameters_array_list:
    def __init__(self, size):
        self.rho_k = np.zeros(size)    
        self.d_1 = np.zeros(size)
        self.d_next = np.zeros(size)
        self.alpha_1eef = np.zeros(size)
        self.Z = np.zeros(size)
        self.Z_new = np.zeros(size)
        # self.Z_ratio = np.zeros(size)
        self.q_t = np.zeros(size)
        self.Delta_H = np.zeros(size)
        self.H_ave = np.zeros(size)
        self.l_11 = np.zeros(size)
        self.d_k = np.zeros(size)
        self.theta = np.zeros(size)
        self.isOK = np.ones(size, dtype=int)
        self.is_d_k_OK = np.ones(size, dtype=int)
        self.__last_it = 0
        
    def set_row(self, isOK, is_d_k_OK, rho_k=-1, d_1=-1, d_next=-1, alpha_1eef=-1, Z=-1, Z_new=-1, q_t=-1, Delta_H=-1, H_ave=-1, l_11=-1, d_k=-1, theta=-1):
        i = self.__last_it
        
        if not(isOK): self.isOK[i] = False
            
        if not(is_d_k_OK):self.is_d_k_OK[i] = False
                
        self.rho_k[i] = rho_k
        self.d_1[i] = d_1
        self.d_next[i] = d_next 
        self.alpha_1eef[i] = alpha_1eef
        self.Z[i] = Z
        self.Z_new[i] = Z_new
        self.q_t[i] = q_t
        self.Delta_H[i] = Delta_H
        self.H_ave[i] = H_ave
        self.l_11[i] = l_11 
        self.d_k[i] = d_k
        self.theta[i] = theta
        self.__last_it += 1

    def set_row_by_split(self, split, isOK, isd_kOK):
        i = self.__last_it
        
        
        self.isOK[i] = isOK
        self.is_d_k_OK[i] = isd_kOK
        
        if split == -1:
            self.rho_k[i] = -1
            self.d_1[i] = -1
            self.alpha_1eef[i] = -1
            self.d_next[i] = -1
            self.Z[i] = -1
            self.Z_new[i] = -1
            self.q_t[i] = -1
            self.Delta_H[i] = -1
            self.H_ave[i] = -1
            self.l_11[i] = -1
            self.d_k[i] = -1
            self.theta[i] = -1
            
            self.__last_it += 1
            return
        
        self.rho_k[i] = split.get_rho_k()
        self.d_1[i] = split.get_d_1()
        self.alpha_1eef[i] = split.get_alpha_1eef()
        
        if isOK == 1:
            self.d_next[i] = split.get_next_d() 
            self.Z[i] = split.get_Z()
            self.Z_new[i] = split.get_Z_new()
            self.q_t[i] = split.get_q_t()
            self.Delta_H[i] = split.get_Delta_H()
            self.H_ave[i] = split.get_H_0ave()
            self.l_11[i] = split.get_l_11() 
            self.d_k[i] = split.get_d_k()
            self.theta[i] = split.get_last_theta()
        else:
            self.d_next[i] = -1
            self.Z[i] = -1
            self.Z_new[i] = -1
            self.q_t[i] = -1
            self.Delta_H[i] = -1
            self.H_ave[i] = -1
            self.l_11[i] = -1
            self.d_k[i] = -1
            self.theta[i] = -1
        
        self.__last_it += 1

class HPC_SPLIT:
    def __init__(
        self,                       # Давление свежего пара                                |
        fp_0,                       # Давление свежего пара после гидр. сопротивления      |
        fh_0,                       # Температура свежего пара                             |
        n,                          # Номинальная частота вращения                         |
        G_0,                        # Расход пара через группу ступеней                    |
        p_z,                        # Давление за группой ступеней                         |
        Z,                          # Предполагаемое число ступеней                        |
        alpha_1eef,                 # Эффективный угол выхода потока из сопловой решетки   |
        etaHPC_oi,                  # Внутренний КПД ЦВД                                   |
        d_1,                        # Средней диаметр первой ступени                       |
        rho_k,                      # Степень реактивности в корне первой ступени          |
        phi = np.mean([0.93, 0.96]),# Коэффициент скорости сопловой решетки первой ступени |
        mu = np.mean([0.95,0.97]),  # Коэффициент расхода сопловой решетки первой ступени  |
        Delta = Q(0.003, "meter"),  # Перекрыша между лопаткиами первой ступени            |
        ):
        # Если что-то пошло не так
        self.__isOK            = True
        
        
        # Входные величины ---------
        #self.fp_0_      = fp_0_
        self.__fp_0       = fp_0
        self.__fh_0       = fh_0
        self.__n          = n
        self.__G_0        = G_0
        self.__p_z        = p_z
        self.__Z          = Z
        self.__alpha_1eef = alpha_1eef
        self.__d_1        = d_1
        self.__rho_k      = rho_k
        self.__phi        = phi
        self.__mu         = mu
        self.__Delta      = Delta
        self.__etaHPC_oi  = etaHPC_oi
        
        # Переменные необходимые в процессе вычисления
        # 2
        self.__fPoint_0_              = None
        self.__fPoint_0               = None
        self.__fs_0                   = None
        
        
        # 3
        self.__theta                  = None
        self.__iteration_number       = None
        self.__percent_difference     = None
        self.__l_11                   = None
        self.__uDIVu_cf               = None
        self.__fH_01                  = None
        self.__h_2t                   = None
        self.  v_2t                   = None
        
        # 4
        self.l_21                   = None
        
        # 5
        self.d_k                    = None
        
        # 6
        self.Point_zt               = None
        self.h_zt                   = None
        self.fH_0                   = None
        self.H_i                    = None
        self.h_z                    = None
        self.Point_2z               = None
        self.v_2z                   = None
        
        # 7
        self.l_2z                   = None
        
        # 8
        self.d_2z                   = None
        
        # 9
        self.theta_1                = None
        self.theta_z                = None
        self.rho_1                  = None
        self.rho_z                  = None
        self.uDIVc_f_1              = None
        self.uDIVc_f_z              = None
        
        # 10
        self.k_for_d                = None
        self.b_for_d                = None
        self.k_for_l                = None
        self.b_for_l                = None
        self.d_vec                  = None
        self.l_vec                  = None
        self.theta_vec              = None
        self.rho_vec                = None
        self.uDIVc_f_vec            = None
        self.K_vec                  = None
        self.H_vec                  = None
        self.H_0ave                 = None
        self.q_t                    = None
        self.Z_new                  = None

        # 11
        self.Delta_H                = None
        
        # 12
        self.H_new_vec              = None
        
        
        # Основной расчет
        # 1 ~
        # 2
        try:
            self.update_full_point_of_income_steam()
        # 3
            self.update_l_11_and_theta_11()
        # 4
            self.update_l_21()
        # 5
            self.update_root_diameter()
        # 6
            self.update_steam_parameters_after_last_stage()
        # 7
            self.update_l_2z()
        # 8
            self.update_d_2z()
        # 9
            self.update_parameters_of_first_and_last_stage()
        # 10
            self.update_heat_transfer()
        # 11
            self.update_Delta_of_heat_transfer()
        # 12
            self.update_heat_transfer_via_Delta()
        except Exception:
            self.__isOK = False
        
        
    # 2. Энтропия и температура пара перед первой ступенью
    def update_full_point_of_income_steam(self):
        #self.fPoint_0_ = TCPv2.ThPoint(p=self.fp_0_, t=self.ft_0)
        self.fPoint_0 = TCPv2.ThPoint(p=self.fp_0, h=self.fh_0)
        self.fs_0     = self.fPoint_0.s()
        self.fh_0     = self.fPoint_0.h()
        
    # 3. Определяем высоту первой нерегулируемой ступени, для этого задаемся величиной обратной веерности
    # theta = d_1/l_1
    def update_l_11_and_theta_11(self):
        
        theta               = 20
        iterations_number   = 0
        percent_difference  = None
        l_11                = None
        uDIVu_cf            = None
        fH_01               = None
        h_2t                = None
        v_2t                = None
        
        while(True):
            # 3.1. Величина обратной веерности
            # ...
    
            # 3.2. Определяем степень реактивности на среднем диаметре
            rho = self.rho_k + 1.8/(theta+1.8)
    
            # 3.3. Определяем оптимальное значение u/c_f
            uDIVu_cf = self.phi*np.cos(self.alpha_1eef)/(2*np.sqrt(1-rho))
    
            # 3.4. Определяем располагаемый теплопеперад по параметрам торможения при оптимальном
            # u/c_f для первой ступени 
            fH_01 = 12300 * np.power((self.d_1 * self.n)/(uDIVu_cf * 50),2)
    
            # 3.5. Определяем теоретическое значение энтальпии за первой ступенью
            h_2t = self.fh_0 - fH_01
    
            # 3.6. Определяем удальный объем пара за первой нерегулируемой ступенью при
            # изоэнтропном процессе расширения по свойствам воды и водяного пара
            Point_2t = TCPv2.ThPoint(h=h_2t, s=self.fs_0)
            v_2t = Point_2t.v()
    
            # 3.7. Определеяем высоту первой нерегулируемой ступени
            l_11 = (self.G_0 * v_2t * uDIVu_cf)/(np.power(np.pi*self.d_1,2) * self.n * np.sqrt(1-rho) * np.sin(self.alpha_1eef) * self.mu)
    
            # 3.8. Определяем окончательное значение обратной веерности и проверяем его
    
            # Проверка условия 
            # Если получившаяся величина больше, чем заданная
            if (self.d_1/l_11 > theta):
                # То вычитаем отношение из единицы
                percent_difference = (1 - theta/(self.d_1/l_11)) 
            # Если величина меньше, чем заданная
            else:
                # То вычитаем из отношения единицу
                percent_difference = (theta/(self.d_1/l_11) - 1)

            # Если условие выполнилось - выходим из цикла
            if (np.abs(percent_difference) < 0.015):
                self.theta               = theta
                self.iterations_number   = iterations_number
                self.percent_difference  = percent_difference
                self.l_11                = l_11 
                self.uDIVu_cf            = uDIVu_cf
                self.fH_01               = fH_01
                self.h_2t                = h_2t
                self.v_2t                = v_2t
                
                break
            # Иначе добавляем итерацию и меняем приближающее значение на найденное в процессе цикла
            else:
                iterations_number += 1
                theta = (self.d_1/l_11)
                
    # 4. Определяем высоту рабочей лопатки первой нерегулируемой ступени
        def update_l_21(self):
            self.l_21 = self.l_11 + self.Delta
    
    # 4. Определяем высоту рабочей лопатки первой нерегулируемой ступени
    def update_l_21(self):
        self.l_21 = self.l_11 + self.Delta
    
    # 5. Определяем корневой диаметр ступени
    def update_root_diameter(self):
        self.d_k = self.d_1 - self.l_21
        
    # 6. Определяем параметры пара за последней ступенью ЦВД
    def update_steam_parameters_after_last_stage(self):
        # 6.1. Значение энтальпии пара при изоэнтропном расширении пара в ЦВД:
        # Термодинамическая точка zt
        self.Point_zt = TCPv2.ThPoint(p=self.p_z, s=self.fs_0)
        # Энтальпия
        self.h_zt = self.Point_zt.h()

        # 6.2. Теоретический перепад на отсек нерегулируемых ступеней ЦВД:
        self.fH_0 = self.fh_0 - self.h_zt

        # 6.3. Действительный теплоперепад на отсек нерегулируемых ступеней ЦВД
        self.H_i = self.fH_0 * self.etaHPC_oi

        # 6.4. Действительное значение энтальпии за ЦВД (за последней ступенью)
        self.h_z = self.fh_0 - self.H_i

        # 6.5 Действительный объем за ЦВД (за последней ступенью)
        # Термодинамическая точка 2z
        self.Point_2z = TCPv2.ThPoint(p=self.p_z, h=self.h_z)
        self.v_2z = self.Point_2z.v()
    
    # 7. Определяем высоту рабочей лопатки последней ступени исходя из
    # из предположения о том, что удельный объем в цилиндре имзеняется 
    # линейно, а так же учитывая закон постоянства корневого диаметра в
    # проточной части ЦВД. Для этого решается квадратное уравнение
    # относительно неизвестной величины
    # l_2z**2 + l_2z*d_k = l_21 * d_21 * v_2z / v_2t
    def update_l_2z(self):
        l_2z = snj_solvers.solve__(
            [
                ('l_2z**2 + l_2z * d_k', 'l_21 * d_1 * v_2z / v_2t'),
            ],
    
            [
                ('d_k',  self.d_k),
                ('l_21', self.l_21),
                ('d_1',  self.d_1),
                ('v_2z', self.v_2z),
                ('v_2t', self.v_2t)
            ],
    
            [
                'l_2z'
            ],
            list_mode=True
        )
        self.l_2z = snj_solvers.return_positive(l_2z[0], l_2z[1])
    
    # 8. Определеяем средний диаметр последней ступени ЦВД
    def update_d_2z(self):
        self.d_2z = self.d_k + self.l_2z
        
    # 9. Определяем основные параметры первой и последней ступени
    def update_parameters_of_first_and_last_stage(self):
        # 9.1. Обратная веерность
        self.theta_1 = (self.l_21 + self.d_k)/self.l_21
        self.theta_z = (self.l_2z + self.d_k)/self.l_2z
        # 9.2. Степень реактивности
        self.rho_1 = self.rho_k + 1.8/(self.theta_1 + 1.8)
        self.rho_z = self.rho_k + 1.8/(self.theta_z + 1.8)
        # 9.3. Оптимальное значение u/c_f
        self.uDIVc_f_1 = self.phi*np.cos(self.alpha_1eef)/(2*np.sqrt(1-self.rho_1))
        self.uDIVc_f_z = self.phi*np.cos(self.alpha_1eef)/(2*np.sqrt(1-self.rho_z))
        
        
    # 10. Производим разбивку теплоперепадов
    class linear_fun_by_number:
        def __init__(self, in_k,in_b):
            self._k = in_k
            self._b = in_b
    
        # получить значение линейной функции
        def get(self, in_z):
            return self._k * in_z + self._b
    
        # получить значение линейной функции без размерности
        def get_no_dim(self, in_z):
            return self.get(in_z).m

        # получить значение линейной функции приведенной к заданной размерности, но без нее
        def get_to_dim_no_dim(self, in_z, dim):
            return self.get(in_z).to(dim).m
    
    def update_heat_transfer(self):    
        # 10.1. Предполагаем, что средний диаметр ступеней и высота лопаток
        # высота лопаток изменяются вдоль ЦВД линейно. Тогда можно построить
        # диаграммы их изменения для цилиндра с количеством ступеней (Z).
        # Таким образом из графиков или по формуле, для каждой ступени имеем
        # величину среднего диаметра и высоты рабочей лопатки.
        self.k_for_d,self.b_for_d = snj_solvers.solve__(
        [
            ('d_1',  'k*z_1 + b'),
            ('d_2z', 'k*z_z + b'),
        ],
        [
            ('d_1', self.d_1),
            ('d_2z', self.d_2z),
            ('z_1', 1),
            ('z_z', self.Z)
        ],
        [
            'k','b'
        ],
        list_mode=True)
        
        self.k_for_l,self.b_for_l = snj_solvers.solve__(
        [
            ('l_21',  'k*z_1 + b'),
            ('l_2z', 'k*z_z + b'),
        ],
        [
            ('l_21', self.l_21),
            ('l_2z', self.l_2z),
            ('z_1', 1),
            ('z_z', self.Z)
        ],
        [
            'k','b'
        ],
        list_mode=True)
        
        # вектора для хранения данных
        self.d_vec = np.array([])
        self.l_vec = np.array([])
        
        # Функция для определения диаметра от номера ступени
        d_fun_handler = self.linear_fun_by_number(self.k_for_d,self.b_for_d)

        # Функция для определения высоты лопатки от номера ступени
        l_fun_handler = self.linear_fun_by_number(self.k_for_l,self.b_for_l)
        
        for i in range(1, self.Z+1):
            self.d_vec = np.append(self.d_vec, d_fun_handler.get(i))
            self.l_vec = np.append(self.l_vec, l_fun_handler.get(i))
            
        # 10.2. Для каждой ступени определяем обратную веерность 
        self.theta_vec = (self.l_vec + self.d_k )/ self.l_vec
        
        # 10.3. Для каждой ступени определяем степень реактивности
        self.rho_vec = self.rho_k + 1.8/(self.theta_vec + 1.8)
        
        # 10.4. Для каждой ступени определяем величину u/c_f
        self.uDIVc_f_vec = self.phi*np.cos(self.alpha_1eef)/(2*np.sqrt(1-self.rho_vec))
        
        # 10.5. Для каждой ступени определяем теплоперепад по статическим параметрам

        # вектор коэффициентов
        self.K_vec = np.empty(self.Z)
        self.K_vec.fill(0.95)
        self.K_vec[0] = 1.0

        self.H_vec = 12300 * (self.d_vec/self.uDIVc_f_vec)**2 * (self.n/50)**2 * self.K_vec
        self.H_vec.ito('kJ/kg')
        
        # 10.6. Определяем среднее значение теплоперепадов
        self.H_0ave = np.mean(self.H_vec)
        
        # 10.7 Определяем коэффициент возврата теплоты:
        q_t_koef = 4.8*10**(-4) * un('kg/kJ')
        self.q_t = q_t_koef*(1 - self.etaHPC_oi)*self.H_0ave*(self.Z-1)/self.Z
        
        # 10.8. Новое значение количества ступеней ЦВД
        self.Z_new = self.fH_0 * (1+self.q_t)/self.H_0ave
        self.Z_new.ito('')

    # 11. Определяем невязку после разбивки теплоперепадов
    def update_Delta_of_heat_transfer(self):
        self.Delta_H = (self.fH_0*(1+self.q_t) - np.sum(self.H_vec))/self.Z
        
    # 12. Уточняем значение теплоперепадов на каждую ступень с учетом невязок
    def update_heat_transfer_via_Delta(self):
        self.H_new_vec = self.H_vec + self.Delta_H
        
def get_auto_number_of_stages_split(
                                    # Давление свежего пара                                |
        fp_0,                       # Давление свежего пара после гидр. сопротивления      |
        fh_0,                       # Температура свежего пара                             |
        n,                          # Номинальная частота вращения                         |
        G_0,                        # Расход пара через группу ступеней                    |
        p_z,                        # Давление за группой ступеней                         |
        alpha_1eef,                 # Эффективный угол выхода потока из сопловой решетки   |
        etaHPC_oi,                  # Внутренний КПД ЦВД                                   |
        d_1,                        # Средней диаметр первой ступени                       |
        rho_k,                      # Степень реактивности в корне первой ступени          |
        phi = np.mean([0.93, 0.96]),# Коэффициент скорости сопловой решетки первой ступени |
        mu = np.mean([0.95,0.97]),  # Коэффициент расхода сопловой решетки первой ступени  |
        Delta = Q(0.003, "meter"),  # Перекрыша между лопаткиами первой ступени            |
        
        #delta_z = 0.1
):
    local_z = 5
    while (True):
        it_split = HPC_SPLIT(
                                  # Давление свежего пара                                |
            fp_0=fp_0,                       # Давление свежего пара после гидр. сопротивления      |
            fh_0=fh_0,                       # Температура свежего пара                             |
            n=n,                          # Номинальная частота вращения                         |
            G_0=G_0,                        # Расход пара через группу ступеней                    |
            p_z=p_z,                        # Давление за группой ступеней                         |
            alpha_1eef=alpha_1eef,                 # Эффективный угол выхода потока из сопловой решетки   |
            etaHPC_oi=etaHPC_oi,                  # Внутренний КПД ЦВД                                   |
            d_1=d_1,                        # Средней диаметр первой ступени                       |
            rho_k=rho_k,                    # Степень реактивности в корне первой ступени          |
            phi = phi,# Коэффициент скорости сопловой решетки первой ступени |
            mu = mu,  # Коэффициент расхода сопловой решетки первой ступени  |
            Delta = Delta,  # Перекрыша между лопаткиами первой ступени            |
            
            Z = local_z
        )
        if not(it_split.__isOK):
            return it_split
        
        #if (it_split.Z_new - local_z < 0):
        #    local_z += 1
        #else:
        #    local_z -= 1
        return it_split
        if not(np.abs(it_split.Z_new-local_z) > 1): 
            #local_z = int(it_split.Z_new.m)
        #else:
            return it_split

# Входные параметры, которые можно варировать
#     d_1         = Q(0.6, "meter"),
#     rho_k       = 0.07
#     alpha_1eef  = Q(9, "deg"),
def calculate_hpc_split_in_ranges(
    d_1,          #(d_1_min, d_1_max, split_value)
    rho_k_range,        #(rho_k_min, rho_k_max, split_value)
    alpha_1eef_range,   #(alpha_1eef_min, alpha_1eef.max, split_value)
    
          #= Q(9, "MPa"),
    fp_0,       #= Q(8.73, "MPa"),
    fh_0,        #= Q(571, "degC"),
    n,           #= Q(90, "1/s"),
    G_0,         #= Q(38.88, "kg/s"),
    p_z,         #= Q(1.8, "MPa"),
    
    #Z           = 2,
    #alpha_1eef  = Q(9, "deg"),
    etaHPC_oi,      #= 0.846,
    phi         = np.mean([0.93, 0.96]),
    mu          = np.mean([0.95,0.97]),
    Delta       = Q(0.003, "meter"),
    
    #d_1         = Q(0.6, "meter"),
    #rho_k       = 0.07,
    #delta_z     = 0.3
):
    '''
    1-dim array
    main_data = np.array(
        []
    )
    
    2-dim array
    main_data = np.array(
        [
            [
                ],
                [    
            ]
        ]
    )
    
    3-dim array 3x3x3 = 27
    main_data = np.array([
        [   # layout 1
            [
                1, 2, 3
            ],
            [
                1, 2, 3
            ],
            [
                1, 2, 3
            ],
        ],  
        [   # layout 2
            [
                1, 2, 3
            ],
            [
                1, 2, 3 
            ],
            [
                1, 2, 3 
            ],
        ], 
        [   # layout 3
            [
                1, 2, 3 
            ],
            [
                1, 2, 3    
            ],
            [
                1, 2, 3  
            ],
        ]  
    ])
            (0,     1,          2)
    np.zeros(d_1,   alpha_1eef, Z)
    
    '''
    # Разбивка по значениям
    #d_1_split           = np.arange(d_1_range[0],        d_1_range[1],          d_1_range[2])
    #d_1_split           = Q(d_1_split, "meter")
    rho_k_split         = np.arange(rho_k_range[0],      rho_k_range[1],        rho_k_range[2])
    rho_k_split         = Q(rho_k_split, "")
    alpha_1eef_split    = np.arange(alpha_1eef_range[0], alpha_1eef_range[1],   alpha_1eef_range[2])
    alpha_1eef_split    = Q(alpha_1eef_split, "deg")
    
    # Количество значений входных переменных
    #d_1_split_number            = len(d_1_split)
    rho_k_split_number          = len(rho_k_split)
    alpha_1eef_split_number     = len(alpha_1eef_split)
    
    #print(d_1_split)
    #print(rho_k_split)
    #print(alpha_1eef_split)
    
    Z_vec = np.zeros((
        rho_k_split_number,
        alpha_1eef_split_number,
        
    ))
    
    l_11_vec = np.zeros((
        rho_k_split_number,
        alpha_1eef_split_number,
    ))
    
    rho_k_vec= np.zeros((
        rho_k_split_number,
        alpha_1eef_split_number,
        #d_1_split_number*alpha_1eef_split_number # Z
    ))
    alpha_1eef_vec= np.zeros((
        rho_k_split_number,
        alpha_1eef_split_number,
        #d_1_split_number*alpha_1eef_split_number # Z
    ))
    
    for rho_k_it in range(0, rho_k_split_number):
        for alpha_1eef_it in range(0, alpha_1eef_split_number):
            loc_split = get_auto_number_of_stages_split(
                fp_0=fp_0,
                
                fh_0=fh_0,
                n=n,
                G_0=G_0,
                p_z=p_z,
                etaHPC_oi=etaHPC_oi,
                rho_k=rho_k_split[rho_k_it],
                alpha_1eef=alpha_1eef_split[alpha_1eef_it],
                d_1=d_1
            )
            
            # loc_split = HPC_SPLIT(
            #     fp_0=fp_0,
            #     fp_0_=fp_0_,
            #     ft_0=ft_0,
            #     n=n,
            #     G_0=G_0,
            #     p_z=p_z,
            #     etaHPC_oi=etaHPC_oi,
            #     Z = 7,
            #     rho_k=rho_k_split[rho_k_it],
                
            #     alpha_1eef=alpha_1eef_split[alpha_1eef_it],
            #     d_1=d_1
            # )
            # if (loc_split.l_21 > Q(13, "mm")):
            #     Z_vec[rho_k_it][alpha_1eef_it] = str(loc_split.Z_new.m) + "(Y)"
            # else:
            #     Z_vec[rho_k_it][alpha_1eef_it] = str(loc_split.Z_new.m) + "(N)"
            if not(loc_split.__isOK):
                l_11_vec[rho_k_it][alpha_1eef_it] = -1.0
                Z_vec[rho_k_it][alpha_1eef_it] = loc_split.Z
            else:
                l_11_vec[rho_k_it][alpha_1eef_it] = loc_split.l_11.m
                Z_vec[rho_k_it][alpha_1eef_it] = loc_split.Z_new.m
            rho_k_vec[rho_k_it][alpha_1eef_it] = rho_k_split[rho_k_it].m
            alpha_1eef_vec[rho_k_it][alpha_1eef_it] = alpha_1eef_split[alpha_1eef_it].m
            
    
    return [rho_k_vec, alpha_1eef_vec, Z_vec, l_11_vec]
    #print('hey')
    #print(number_of_d_1)



from io import open
def save_split_data_in_tex(d_1,  rho_k_vec, alpha_1eef_vec, Z_vec):
    
    def swap_br(in_str:str):
        return in_str.replace('[', "{").replace(']', "}")
    
    rho_k_vec_table = []
    for i in range(0,len(rho_k_vec)):
        rho_k_vec_table.append(round(rho_k_vec[i][0], 5))
    #print("rho_k_vec_table:", rho_k_vec_table)
    
    alpha_1eef_vec_table = []
    for i in range(0,len(alpha_1eef_vec[0])):
        alpha_1eef_vec_table.append(round(alpha_1eef_vec[0][i], 5))
    #print("alpha_1eef_vec_table:", alpha_1eef_vec_table)
    
    with open("tex/{}.tex".format("d_1_it1_1"), 'w',encoding="UTF-8") as file:
  
  
        file.write(r"\documentclass[12pt]{article}"+"\n")
        file.write(r"\usepackage[utf8]{inputenc}"+"\n")
        file.write(r"\usepackage{mathtext}" + '\n')
        file.write(r"\usepackage[T2A]{fontenc}"+"\n")
        file.write(r"\usepackage[english, russian]{babel}" +"\n")
        
        file.write(r"\begin{document}"+"\n")
        #file.write(r"Данный файл создан для разбавки с диаметром d={}\\".format(d_1)+"\n")
        
        number_of_columns = "c|"
        for i in range(0, len(rho_k_vec_table)):
            number_of_columns += 'c'
        
    #     '''
    #     \hline
    #         \multicolumn{4}{c}{$d_{1}=10\  м$} \\ 
    #     \hline
    #         $\alpha_{1эф}$ & \multicolumn{3}{|c}{$\rho$} \\ 
    #     \hline
    #         9 & 1 & 1 & 1 \\
    #         10 & 1 & 1 & 1 \\
    #         11 & 1 & 1 & 1 \\
    #         12 & 1 & 1 & 1 \\
    #         13 & 1 & 1 & 1 \\
    #     '''
        
        file.write(swap_br(r"\begin[tabular][{}]".format(number_of_columns)) + "\n")
        file.write(r"\hline" + '\n')
        file.write('\t' + swap_br(r'\multicolumn[{}][{}][{}] \\'.format(
            str(len(rho_k_vec_table) + 1),
            'c',
            swap_br('$d_[1]={} \ м$'.format(d_1))
        )) + '\n')
        file.write(r"\hline" + '\n')
        
        # & \multicolumn{3}{|c}{$\rho$} \\ 
        file.write('\t' + swap_br(r' & \multicolumn[{}][{}][{}] \\'.format(
            str(len(rho_k_vec_table)),
            '|c',
            r"$\rho$"
        )) + '\n')
        
    #     # \cline{2-4}
        file.write('\t' +  swap_br(r"\cline[2-{}]".format(
            str(len(rho_k_vec_table) + 1)
        )) + '\n')
        
        #\raisebox{1.5ex}[0cm][0cm]{$\alpha_{1эф}$} & 1 & 2 & 3 \\ 
        rho_arange_str = r" & "
        for i in range(0, len(rho_k_vec_table)):
            rho_arange_str += str(rho_k_vec_table[i])
            if not(i == len(rho_k_vec_table) - 1): rho_arange_str += " & "
        rho_arange_str += r"\\"    
        
        file.write('\t' + swap_br(r"\raisebox[1.5ex]") + r"[0cm][0cm]" + r"{$\al"+ "pha_{1эф}$}" + rho_arange_str + '\n')
        
        file.write(r"\hline" + '\n')
        
        

        
        for i in range(0, len(alpha_1eef_vec_table)):
            data_row_str = r"" + str(alpha_1eef_vec_table[i]) + r" & " + '\t'
            for j in range(0, len(rho_k_vec_table)):
                # \textbf
                if ((abs(int(Z_vec[j][i]) - Z_vec[j][i])) < 0.1):
                    data_row_str += r"\textbf{"
                    data_row_str += str(round(Z_vec[j][i], 4))
                    data_row_str += r"}"
                else:
                    data_row_str += str(round(Z_vec[j][i], 4))
                if not(j == (len(rho_k_vec_table) - 1)): data_row_str += r" & " + '\t'
            data_row_str += r"\\"
            file.write('\t' + data_row_str + '\n')
        # '''
        # \hline
        # \multicolumn{4}{c}{$d_{1}=10\  м$} \\ 
        # \hline
        #     & \multicolumn{3}{|c}{$\rho$} \\ 
        #     \cline{2-4}
        #     \raisebox{1.5ex}[0cm][0cm]{$\alpha_{1эф}$} & 1 & 2 & 3 \\ 
        # \hline
        #     9 & 1 & 1 & 1 \\
        #     10 & 1 & 1 & 1 \\
        #     11 & 1 & 1 & 1 \\
        #     12 & 1 & 1 & 1 \\
        #     13 & 1 & 1 & 1 \\'''
        #file.write(r"\multicolumn[{}][{}][{}]\\".format(len(str(len(rho_k_vec)-1))))
        
        #file.write(r"$\rho$ \\ " + "\n")
        #file.write(r"MMMMMMMM \= MMMMMMMM \= MMMMMMMM \kill" + "\n")
        #file.write(r"1 \>2 \> 3" + "\n")
        #file.write(r""+"\n")
        
        
        file.write(r"\end{tabular}" + "\n")
        file.write(r"\end{document}"+"\n")

def generate_simple_tex_file(path):
    if '.tex' in path: path = path.replace('.tex', '')
    with open(path  + ".tex", 'w', encoding="UTF-8") as file:
        file.write(r"\documentclass[20pt]{article}"+"\n")
        file.write(r"\usepackage[left=3cm,right=1cm,top=2cm,bottom=2cm,bindingoffset=0cm]{geometry}")
        file.write(r"\usepackage[utf8]{inputenc}"+"\n")
        file.write(r"\usepackage{mathtext}" + '\n')
        file.write(r"\usepackage[T2A]{fontenc}"+"\n")
        file.write(r"\usepackage[english, russian]{babel}" +"\n")
        file.write(r"\begin{document}"+"\n")
        file.write(r"\begin{center}" + "\n")
        file.write(r"%!LOCAL_BEGIN!" + "\n")
        
        file.write(r"\end{center}" + '\n')
        file.write(r"\end{document}"+"\n")

def add_smt_to_existing_file(path, data:str):
    lines = None
    with open(path, 'r', encoding="UTF-8") as file:
        lines = file.readlines()
    with open(path, 'w', encoding="UTF-8") as file:
        for i in lines:
            if ("%!LOCAL_BEGIN!" in i):
                file.write(data)
                file.write("\n" + "%!LOCAL_BEGIN!" + "\n")
            else:
                file.write(i)
        
        # cursor_pos = 0
        # for item in file:
        #     if ("!LOCAL_BEGIN!" in item):
                
        #         break
        #     else:
        #         cursor_pos += len(item) + 1
        # file.seek(cursor_pos)
        # file.write(data)
                
        

def generate_table_of_split_data_in_tex(d_1, rho_k_vec, alpha_1eef_vec, Z_vec, l_11_vec):
    def swap_br(in_str:str):
        return in_str.replace('[', "{").replace(']', "}")
    
    rho_k_vec_table = []
    for i in range(0,len(rho_k_vec)):
        rho_k_vec_table.append(round(rho_k_vec[i][0], 5))
    #print("rho_k_vec_table:", rho_k_vec_table)
    
    alpha_1eef_vec_table = []
    for i in range(0,len(alpha_1eef_vec[0])):
        alpha_1eef_vec_table.append(round(alpha_1eef_vec[0][i], 5))
    #print("alpha_1eef_vec_table:", alpha_1eef_vec_table)
    
    
    main_data_string = r""
    
        
        
    number_of_columns = "c|"
    for i in range(0, len(rho_k_vec_table)):
        number_of_columns += 'c'
        
    main_data_string+=(swap_br(r"\begin[tabular][{}]".format(number_of_columns)) + "\n")
    main_data_string+=(r"\hline" + '\n')
    main_data_string+=('\t' + swap_br(r'\multicolumn[{}][{}][{}] \\'.format(
        str(len(rho_k_vec_table) + 1),
        'c',
        swap_br('$d_[1]={} \ м$'.format(d_1))
    )) + '\n')
    main_data_string+=(r"\hline" + '\n')
        
    # & \multicolumn{3}{|c}{$\rho$} \\ 
    main_data_string+=('\t' + swap_br(r' & \multicolumn[{}][{}][{}] \\'.format(
        str(len(rho_k_vec_table)),
        '|c',
        r"$\rho$"
    )) + '\n')
        
    #     # \cline{2-4}
    main_data_string+=('\t' +  swap_br(r"\cline[2-{}]".format(
        str(len(rho_k_vec_table) + 1)
    )) + '\n')
        
        #\raisebox{1.5ex}[0cm][0cm]{$\alpha_{1эф}$} & 1 & 2 & 3 \\ 
    rho_arange_str = r" & "
    for i in range(0, len(rho_k_vec_table)):
        rho_arange_str += str(rho_k_vec_table[i])
        if not(i == len(rho_k_vec_table) - 1): rho_arange_str += " & "
    rho_arange_str += r"\\"    
        
    main_data_string+=('\t' + swap_br(r"\raisebox[1.5ex]") + r"[0cm][0cm]" + r"{$\al"+ "pha_{1эф}$}" + rho_arange_str + '\n')
        
    main_data_string+=(r"\hline" + '\n')
        
    for i in range(0, len(alpha_1eef_vec_table)):
        data_row_str = r"" + str(alpha_1eef_vec_table[i]) + r" & " + '\t'
        for j in range(0, len(rho_k_vec_table)):
            # \textbf
            if (l_11_vec[j][i] == -1.0):
                data_row_str += r"Нет"
            else:
            
                if (((abs(int(Z_vec[j][i]) - Z_vec[j][i])) < 0.1)and(l_11_vec[j][i] >= 0.013)):
                    data_row_str += r"\textbf{"
                    data_row_str += str(round(Z_vec[j][i], 2)) + r"|" + str(round(l_11_vec[j][i], 4))
                    data_row_str += r"}"
                else:
                    data_row_str += str(round(Z_vec[j][i], 2)) + r"|" + str(round(l_11_vec[j][i], 4))
            if not(j == (len(rho_k_vec_table) - 1)): data_row_str += r" & " + '\t'
        
        data_row_str += r"\\"
        main_data_string+=('\t' + data_row_str + '\n')
    main_data_string+=(r"\end{tabular}" + "\n")
    return main_data_string


# fff = HPC_SPLIT(
#     d_1 = Q(0.6,'m'),        
#     rho_k=Q(0.05, ""),
#     alpha_1eef_range=Q(9, 'deg'),       
#     fp_0        = Q(8.73, "MPa"),
#     fh_0        = Q(3563.4, "kJ/kg"),
#     n           = Q(95, "1/s"),
#     G_0         = Q(38.88, "kg/s"),
#     p_z         = Q(1.8, "MPa"),
#     etaHPC_oi   = 0.846,
#     mu = Q(0.98, ''),
#     phi = Q(0.98, ''),
#     Z = 5
# )

# ddd = calculate_hpc_split_in_ranges(
#         d_1 = Q(0.6, "meter"),        
#         rho_k_range=(0.03,0.10, 0.01),      
#         alpha_1eef_range=(9,15, 1),       
#         fp_0_       = Q(9, "MPa"),
#         fp_0        = Q(8.73, "MPa"),
#         ft_0        = Q(571, "degC"),
#         n           = Q(90, "1/s"),
#         G_0         = Q(38.88, "kg/s"),
#         p_z         = Q(1.8, "MPa"),
#         etaHPC_oi   = 0.846
#     )

# d_split = np.arange(1, 1.2, 0.1)
# d_split = Q(d_split, "meter")
# splits_data = []

# for i in range(0, len(d_split)):
#     splits_data.append(calculate_hpc_split_in_ranges(
#         d_1 = d_split[i],        
#         rho_k_range=(0.05,0.07,0.01),      
#         alpha_1eef_range=(12,20, 1),       
        
#         fp_0        = Q(0.582, "MPa"),
#         fh_0        = Q(2860.1, "kJ/kg"),
        
#         n           = Q(50, "1/s"),
#         G_0         = Q(84.96, "kg/s"),
#         p_z         = Q(0.16, "MPa"),
#         etaHPC_oi   = 0.892
#     ))



# a = []
# d_split = np.arange(0.6, 0.9, 0.1)
# d_split = Q(d_split, "meter")
# for i in range(0, len(d_split)):
#     a.append(calculate_hpc_split_in_ranges(
#         d_1 = d_split[i],        
#         rho_k_range=(0.05,0.07, 0.01),      
#         alpha_1eef_range=(9,20, 1),       
#         fp_0        = Q(0.485, "MPa"),
#         fh_0        = Q(2867, "kJ/kg"),
#         n           = Q(90, "1/s"),
#         G_0         = Q(47.93, "kg/s"),
#         p_z         = Q(0.25, "MPa"),
#         etaHPC_oi   = 0.887
#     ))

# Point_0_ = TCPv2.ThPoint(p=Q(9, "MPa"),t=Q(571, "degC"))
# Point_0  = TCPv2.ThPoint(p=Q(8.73, "MPa"), h=Point_0_.h())

# eee = calculate_hpc_split_in_ranges(
#         d_1 = 0.6*un('meter'),        
#         rho_k_range=(0.03,0.07, 0.01),      
#         alpha_1eef_range=(9,20, 1),       
#         fp_0        = Q(8.73, "MPa"),
#         fh_0        = Point_0.h(),
#         n           = Q(90, "1/s"),
#         G_0         = Q(38.88, "kg/s"),
#         p_z         = Q(1.8, "MPa"),
#         etaHPC_oi   = 0.846
#     )
# # save_split_data_in_tex(
# #     0.6,
# #     ddd[0],
# #     ddd[1],
# #     ddd[2]
# # )
   
# ddd1 = (generate_table_of_split_data_in_tex(
#     0.6,
#     ddd[0],
#     ddd[1],
#     ddd[2]
# ))

# eee1 = (generate_table_of_split_data_in_tex(
#     0.7,
#     eee[0],
#     eee[1],
#     eee[2]
# ))

# generate_simple_tex_file('tex/cccc.tex')
# add_smt_to_existing_file('tex/cccc.tex', ddd1)
# add_smt_to_existing_file('tex/cccc.tex', eee1)
#     #print(main_data[0][0][0])




# calculate_hpc_split_in_ranges(
#     d_1_range=(0.7,0.9, 0.1),            #(d_1_min, d_1_max)
#     rho_k_range=(0.03,0.07, 0.1),        #(rho_k_min, rho_k_max)
#     alpha_1eef_range=(9,14, 1),        #(alpha_1eef_min, alpha_1eef.max)
    
#     fp_0_       = Q(9, "MPa"),
#     fp_0        = Q(8.73, "MPa"),
#     ft_0        = Q(571, "degC"),
#     n           = Q(90, "1/s"),
#     G_0         = Q(38.88, "kg/s"),
#     p_z         = Q(1.8, "MPa"),
    
#     #Z           = 2,
#     #alpha_1eef  = Q(9, "deg"),
#     etaHPC_oi   = 0.846,
#     #phi         = np.mean([0.93, 0.96]),
#     #mu          = np.mean([0.95,0.97]),
#     #Delta       = Q(0.003, "meter"),
    
#     #d_1         = Q(0.6, "meter"),
#     #rho_k       = 0.07,
#     #delta_z     = 0.3
# )
# test_split = HPC_SPLIT(
#     fp_0_       = Q(9, "MPa"),
#     fp_0        = Q(8.73, "MPa"),
#     ft_0        = Q(571, "degC"),
#     n           = Q(90, "1/s"),
#     G_0         = Q(38.88, "kg/s"),
#     p_z         = Q(1.8, "MPa"),
#     Z           = 9,
#     alpha_1eef  = Q(9, "deg"),
#     etaHPC_oi   = 0.846,
#     d_1         = Q(0.6, "meter"),
#     rho_k       = 0.07
# )