#   NGV - nozzle guide vanes
#   RB - rotor blades

'''
    Необходимо определить частные характеристики:
    
    1. rel_t - относительный оптимальный шаг решетки
    2. inst_angle - установочный угол решетки
    3. (возможно) b_model - модельную хорду
    4. a - горло
    5. Delta_in - радиус входной кромки
    6. Delta_out - радиус выходной кромки
'''

import numpy as np
from . import ThPoint

def calculate_lambda(inspeed:float, inThPoint:ThPoint) -> float:
    return inspeed / (np.sqrt(inThPoint.k(0) * inThPoint.p(0) * 10**6 * inThPoint.v(0)))

# rel_t
def calculate_rel_t(in_angle:float, out_angle:float, end_speed:float, end_ThPoint:ThPoint) -> float:
    K = np.sin(np.deg2rad(in_angle)) / np.sin(np.deg2rad(out_angle))
    theta = 180 - (in_angle + out_angle)
    
    if (K < 1.5):
        t_rel_opt_0 = ((1.727 / K) - 0.869) * (1 / np.cbrt(theta)) - (1.71 / K) + 1.604
    else:
        t_rel_opt_0 = (0.327 / np.power(K, 0.371) / np.cbrt(theta)) - (0.994 / np.power(K, 0.395)) + 1.314

    # Приведенная скорость на выходе из решетки 
    lambda_02 = calculate_lambda(end_speed, end_ThPoint)
    
    # Абсолютная поправка по приведенной скорости
    delta_t_rel_opt = -0.625 * lambda_02 ** 2 + 0.48 * lambda_02 + 0.016
    
    # TODO: Варируемое 
    k_critical = 1.04
    
    return t_rel_opt_0 * k_critical * (1 + delta_t_rel_opt)
    
    
    
    
# inst_angle
def calculate_inst_angle(in_angle: float, out_angle:float) -> float:
    return 57.84 - 0.3928 * in_angle + 0.8221 * out_angle

# a
def calculate_throat(end_speed:float, blade_step:float, end_ThPoint:ThPoint) -> float:
    lambda_02 = calculate_lambda(inspeed=end_speed, inThPoint=end_ThPoint)
    k = end_ThPoint.k(0)
    q = np.power((k + 1) / 2, 1/k - 1) * lambda_02 * np.power(1 - lambda_02 ** 2 * (k - 1) / (k + 1), 1/k - 1)
    
    return blade_step * q


# b_model

# def calculate_c_max(chord):
#     return np.sqrt(2000 * in_ThPoint.h(0))

# delta_in
def calculate_edge_radiuses(chord, type):
    r_2 = calculate_r_2(chord)
    r_1 = calculate_r_1(r_2, type)
    return r_1, r_2

def calculate_r_1(r_2, type:str):
    if type == "r":
        return r_2 * 2
    else:
        return r_2 * 10
    
    
    return 0.0527 * np.sin(np.deg2rad(in_angle)) + 0.0071 * np.sin(np.deg2rad(out_angle)) + 0.237 * c_max + 0.18 * r_2 - 0.053


# delta_out
def calculate_r_2(chord):
    return np.mean([0.01, 0.05]) / 2 * chord