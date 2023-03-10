

from . import snj_solvers
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import pint
from pint import get_application_registry, Quantity
un = get_application_registry()
Q = Quantity

import pandas as pd
import numpy as np

# Термодинамические свойства воды и водяного пара
from . import TCPv2
from itertools import product

def HPC_SPLIT_RANGES(                    
    fp_0,                                   # Давление свежего пара после гидр. сопротивления      |
    fh_0,                                   # Энтальпия свежего пара                               |
    n,                                      # Номинальная частота вращения                         |
    G_0,                                    # Расход пара через группу ступеней                    |
    p_z,                                    # Давление за группой ступеней                         |
    etaHPC_oi,                              # Внутренний КПД ЦВД                                   |
    
    alpha_1eef_range = Q([9, 20, 1], "deg"),   # Эффективный угол выхода потока из сопловой решетки   |
    d_1_range = Q([0.6, 1.4, 0.1], "m"),         # Средней диаметр первой ступени                       |
    rho_k_range = Q([0.03, 0.07, 0.01], ""),       # Степень реактивности в корне первой ступени          |
    #phi = np.mean([0.93, 0.96]),           # Коэффициент скорости сопловой решетки первой ступени |
    #mu = np.mean([0.95,0.97]),             # Коэффициент расхода сопловой решетки первой ступени  |
    #Delta = Q(0.003, "meter"),             # Перекрыша между лопаткиами первой ступени            |
    #Z = Q(6,""),                           # Предполагаемое число ступеней                        |
    MODE = "OVERLOAD"                       # Режим расчета ступеней, перегруженный "OVERLOAD"     |                          # или недогруженный "UNDERLOAD"  
):
    if (isinstance(d_1_range,pint.Quantity)):
        d_1_vec_ = np.arange(*(d_1_range.m))
        #d_1_vec = Q(d_1_vec, "m")
    else:
        d_1_vec_ = np.arange(*(d_1_range))
        #d_1_vec = Q(d_1_vec, "m")
    if (isinstance(alpha_1eef_range, pint.Quantity)):
        alpha_1eef_vec_ = np.arange(*(alpha_1eef_range.m))
        #alpha_1eef_vec = Q(alpha_1eef_vec, "deg")
    else:
        alpha_1eef_vec_ = np.arange(*(alpha_1eef_range))
        #alpha_1eef_vec = Q(alpha_1eef_vec, "deg")
    
    if (isinstance(rho_k_range, pint.Quantity)):
        rho_k_vec_ = np.arange(*(rho_k_range.m))
        #rho_k_vec = Q(rho_k_vec, "")
    else:
        rho_k_vec_ = np.arange(*(rho_k_range))
        #rho_k_vec = Q(rho_k_vec, "")
    
    
      
    Z_vec = Q(np.zeros((
        len(d_1_vec_) *
        len(alpha_1eef_vec_) *
        len(rho_k_vec_)
    )),"")
    Z_new_vec = Q(np.zeros((
        len(d_1_vec_) *
        len(alpha_1eef_vec_) *
        len(rho_k_vec_)
    )),"")
    l_11_vec = Q(np.zeros((
        len(d_1_vec_) *
        len(alpha_1eef_vec_) *
        len(rho_k_vec_)
    )),"m")
    d_1_vec = Q(np.zeros((
        len(d_1_vec_) *
        len(alpha_1eef_vec_) *
        len(rho_k_vec_)
    )),"m")
    alpha_1eef_vec = Q(np.zeros((
        len(d_1_vec_) *
        len(alpha_1eef_vec_) *
        len(rho_k_vec_)
    )),"deg")
    rho_k_vec = Q(np.zeros((
        len(d_1_vec_) *
        len(alpha_1eef_vec_) *
        len(rho_k_vec_)
    )),"")
    next_d_vec = Q(np.zeros((
        len(d_1_vec_) *
        len(alpha_1eef_vec_) *
        len(rho_k_vec_)
    )),"m")
    isOK_vec = np.zeros((
        len(d_1_vec_) *
        len(alpha_1eef_vec_) *
        len(rho_k_vec_)
    ))
    
    
    
    counter = 0
    for it in product(d_1_vec_,alpha_1eef_vec_,rho_k_vec_):
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
            rho_k                   =Q(it[2], '')
        )
        if (local_split.isOK):
            Z_vec[counter] = local_split.get_Z()
            Z_new_vec[counter] = local_split.get_Z_new()
            l_11_vec[counter] = local_split.get_l_11()
            next_d_vec[counter] = local_split.get_next_d()
            
            d_1_vec[counter] = Q(it[0],'m')
            alpha_1eef_vec[counter] = Q(it[1], 'deg')
            rho_k_vec[counter] = Q(it[2], '')
            
            isOK_vec[counter] = True
        else:
            Z_vec[counter] = Q(-1, '')
            Z_new_vec[counter] = Q(-1, '')
            l_11_vec[counter] = Q(-1, 'm')
            next_d_vec[counter] = Q(-1, 'm')
            
            d_1_vec[counter] = Q(it[0],'m')
            alpha_1eef_vec[counter] = Q(it[1], 'deg')
            rho_k_vec[counter] = Q(it[2], '')
            
            isOK_vec[counter] = False
        counter += 1
    
    data = {
        'n':np.full(len(d_1_vec), n),
        'd_1':d_1_vec,
        'alpha_1eef':alpha_1eef_vec,
        'rho_k':rho_k_vec,
        'Z':Z_vec,
        'Z_new':Z_new_vec,
        'l_11':l_11_vec,
        'd_next':next_d_vec,
        'OK?':isOK_vec
    }
    return pd.DataFrame(data)

# columns=['n', 'd_1', 'alpha_1eef', 'rho_k']
    df = pd.DataFrame(data)
    
    return [d_1_vec, alpha_1eef_vec, rho_k_vec, Z_vec, Z_new_vec, l_11_vec, isOK_vec]

class HPC_SPLIT_2:
    def __init__(
        self,                       
        fp_0,                       # Давление свежего пара после гидр. сопротивления      |
        fh_0,                       # Энтальпия свежего пара                               |
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
        Z = Q(6,""),                      # Предполагаемое число ступеней                        |
        MODE = "OVERLOAD"           # Режим расчета ступеней, перегруженный "OVERLOAD"     |
                                    # или недогруженный "UNDERLOAD"                        |
    ):
        self.__fp_0 = fp_0
        self.__fh_0 = fh_0
        self.__n    = n
        self.__G_0  = G_0
        self.__p_z  = p_z
        self.__Z    = Z
        self.__alpha_1eef = alpha_1eef
        self.__etaHPC_oi = etaHPC_oi
        self.__d_1 = d_1
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
        
        # Если в процессе расчета высоты первой лопатки произошло исключение вариант не считается)
        self.isOK                   = True
        
        
        self.__calculate()

        
        
        
        
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
            d_k = self.__d_1 - l_21
            d_k.ito('m')
            
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
            l_2z = (-d_k + np.sqrt(d_k**2 + 4*l_21*self.__d_1*v_2z/v_2t))/2
            
            # 8 Средний диаметр последней ступени группы
            d_2z = d_k + l_2z
            
            # 9.1 Обратная вверность в первой и последней ступени
            theta_1 = (l_21 + d_k)/l_21
            theta_z = (l_2z + d_k)/l_2z
            # 9.2 Степень реактивности в первой и последней ступени
            rho_1 = self.__rho_k + 1.8/(theta_1 + 1.8)
            rho_z = self.__rho_k + 1.8/(theta_z + 1.8)
            # 9.3 Оптимальное значение u/c_ф
            uDIVu_1 = self.__phi*np.cos(self.__alpha_1eef)/(2*np.sqrt(1-rho_1))
            uDIVu_z = self.__phi*np.cos(self.__alpha_1eef)/(2*np.sqrt(1-rho_z))
            
            ## 10.1.1, 10.1.2
            # Коэффициент k для функции диаметра
            k_d = (d_2z - self.__d_1)/(self.__Z.m-1)
            b_d = (self.__d_1*self.__Z.m - d_2z)/(self.__Z.m-1)

            # Коэффициент k для функции длины лопатки
            k_l = (l_2z - l_21)/(self.__Z.m-1)
            b_l = (l_21*self.__Z.m - l_2z)/(self.__Z.m-1)

            # Функция для определения диаметра по номеру ступени
            def d(z_i):
                return k_d*z_i + b_d

            # Функция для определения длины лопатки по номеру ступени
            def l(z_i):
                return k_l*z_i + b_l

            # Перечислим номера ступеней в данном векторе
            self.__stages_number_vec = np.arange(1, self.__Z.m+1, 1)

            # Диаметры каждой ступени
            self.__d_vec = d(self.__stages_number_vec)

            # Длины сопловых лопаток каждой ступени
            self.__l_vec = l(self.__stages_number_vec)
            
            # 10.2
            # Обратная вверность для каждой ступени
            self.__theta_vec = (self.__l_vec + d_k)/self.__l_vec

            # 10.3
            # Степень реактивности на среднем диаметре для каждой ступени
            self.__rho_vec = self.__rho_k + 1.8/(self.__theta_vec + 1.8)
            
            # 10.4. Для каждой ступени определяем величину u/c_f
            uDIVc_f_vec = self.__phi*np.cos(self.__alpha_1eef)/(2*np.sqrt(1-self.__rho_vec))

            # 10.5 Теплоперепад по статическим параметрам для каждой ступени
            # Вектор коэффициентов K_i
            self.__K_vec = np.full(self.__Z.m, 0.95)
            self.__K_vec[0] = 1.0

            self.__H_vec = 12300 * (self.__d_vec/uDIVc_f_vec)**2 * (self.__n/50)**2 * self.__K_vec
            self.__H_vec.ito('kJ/kg')
            
            # 10.6 Среднее значение теплоперепада за группу ступеней
            self.__H_0ave = np.mean(self.__H_vec)
            
            # 10.7 Коэффициент возврата теплоты
            q_t_k = Q(4.8*10**(-4), 'kg/kJ')
            self.__q_t = q_t_k * (1 - self.__etaHPC_oi)*self.__H_0ave * (self.__Z.m-1)/self.__Z.m
            
            # 10.8 Уточненное количество ступеней группы
            self.__Z_new = fH_0/self.__H_0ave * (1+self.__q_t)
            self.__Z_new.ito('')

            # 11 Величина распределения теплоперепадов дробной ступени на остальные
            self.__Delta_H = (fH_0*(1+self.__q_t)/self.__Z.m) - self.__H_0ave
            
            # 12 Уточненные теплоперепады
            self.__H_new_vec = self.__H_vec + self.__Delta_H
            
            # Диаметр первой ступени следующего цилиндра
            self.__next_d = self.__d_vec[len(self.__d_vec)-1] + (self.__d_vec[1]-self.__d_vec[0])
            
            # Для автоматическогов выбора ступеней
            if ((self.__MODE == "OVERLOAD")and not(abs(int(self.__Z_new.m) - self.__Z.m) < 1)):
                self.__Z = Q(int(self.__Z_new), "")
                self.__calculate()
            #elif (((self.__MODE == "OVERLOAD")and (int(self.__Z_new.m) == self.__Z.m + 1))):
                #pass# not(int(self.__Z_new.m)+1 == self.__Z.m)
            if ((self.__MODE == "UNDERLOAD")and not(int(self.__Z_new.m) - self.__Z.m > 1)):
                self.__Z = Q(int(self.__Z_new) + 1, "")
                self.__calculate()
            #else:
                #pass
                
            
        except Exception as s:
            self.isOK = False
            #print(s)
    
    # getters
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
        self.isOK            = True
        
        
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
            self.isOK = False
        
        
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
        if not(it_split.isOK):
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
            if not(loc_split.isOK):
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


'''
\begin{tabbing}
    \hspace{l} \= MMMM \= MMMM \kill
    #1\>abc\>dddd\\
    #2\>abc\>dddd\\
    #3\>abc\>dddd\\
\end{tabbing}
'''

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