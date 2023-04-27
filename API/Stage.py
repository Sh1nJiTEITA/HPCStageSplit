# Используемые библиотеки

# Общие
import numpy as np
import pandas as pd
from pint import Quantity as Q

# Собственные
from . import ThPoint
from . import calculate_losses






class Stage:
    """
        Класс "Ступень" создан для рассчета ступеней паровых турбин.
        
    """
    
    def __init__(
        self, 
        p_0,            
        t_0,
        G_0,            # [кг/с] Массовый расход через ступень
        d_hub, 
        n, 
        reaction, 
        alpha_0,
        alpha_1eef, 
        H_0, 
        c_0,  
        Delta_pr,       # [м] Величина перекрыши
        
        extra_param:dict = {},
        
        kappa_vs = 0.5,
    
        **kwargs
    ):
        ''' 
        Инициализация класса:
        p0          - давление перед ступенью
        t0          - темпертура перед ступенью
        d_hub       - корневой диаметр
        n           - частота вращения
        reaction    - степень реактивности на среднем диаметре
        alpha1      - угол выхода из сопловой решетки
        H0          - теплоперепад по статическим параметрам
        c0          - скорость на входе
        alpha0      - угол входа в ступень
        kappa_vs    - каппа вс
        
        extra_param - словарь с дополнительными параметрами для ступени. 
        Список дополнительных параметров:
        
        "u/c_f": float()      - отношение, если ступень рассчитывается не на оптимальные
        параметры.
        
        
        
        '''     

        self.__is_stage_optimal:bool = True # Флаг оптимальности ступени

        # Дополнительные параметры
        if not len(extra_param) == 0:
            for param in extra_param:
                # Если задано u/c_f, то ступень не оптимальная
                if (param == "u/c_f"):
                    self.__is_stage_optimal = False
                    self.__X = extra_param.get(param)    
        
        # ,,,
        #self.__h_1t         = None      # [кДж/кг]   Теоретическая статическая энтальпия на входе в сопловую решетку
        #self.__h_1          = None      # [кДж/кг]   Действительная статическая энтальпия на входе в сопловую решетку
        #self.__h_2t         = None      # [кДж/кг]   Теоретическая энтальпия на выходе из рабочей решетки
        #self.M_i            = None
        #self.__X        = None          # []         Значение (u/c_f)_отп
        #self.__fp_0         = None      # [МПа]      Давление полного торможения на входе в ступень
        #self.__fh_0         = None      # [кДж/кг]   Энтальпия полного торможения на входе в ступень
        #self.__fs_0         = None      # [кДж/кг/К] Энтропия полного торможения на входе в ступень
        #self.__p_1          = None      # [МПа]     Давление на выходе из сопловой решетки
        #self.__p_2          = None      # [МПа]     Статическое давление на выходе из ступени
        #self.__p_0 = p_0                # [МПа]      Статическое давление на входе в ступень
        #self.__t_0 = t_0                # [degC]     Статическая температура на входи в ступень
        #self.__v_1t         = None      # [м^3/кг]   Теоретический удельный объем после сопловой решетки
       
        
        # Термодинамические точки ----- # ----------------------------------------------------------------------------------
                                           
        self.__Point_0 = ThPoint(p=p_0,t=t_0) # ...  Т.Т. статических параметров на входе в ступень
        self.__fPoint_0     = None      # ...        Т.Т. полных параметров на входе в ступень
        self.__Point_1t     = None      # ...        Т.Т. теоретических статических параметров после сопловой решетки
        self.__Point_1      = None      # ...        Т.Т. действительных параметров после сопловой решетки
        self.__fPoint_1     = None      # ...        T.T. полных параметров перед рабочей решеткой
        self.__Point_2t     = None      # ...        Т.Т. теоретических статических параметров после рабочей решетки
        self.__Point_2      = None      # ...        Т.Т. действительных параметров после рабочей решетки
        self.__Point_2tt    = None      # ...        Т.Т. теоретических статических параметров после ступени (2't)
        self.__Point_outt   = None      # ...        Т.Т. действительных статических параметров на выходе из ступени
        self.__Point_out    = None      # ...        Т.Т. действительных полных параметров на выхходе из ступени
        
        # Коэффициенты скорости ------- # ----------------------------------------------------------------------------------
                                        
        self.__phi          = None      # []         Коэффициент скорости для соплвоой решетки
        self.__psi          = None      # []         Коэффициент скорости для рабочей решетки
                                         
        # Коэффициенты расхода -------- # ----------------------------------------------------------------------------------
                                        
        self.__mu_1         = None      # []         Коэффициент расхода для сопловой решетки
        self.__mu_2         = None      # []         Коэффициент расхода для рабочей решетки
                                        
        # Теплоперепады решеток ------- # ----------------------------------------------------------------------------------
        
        self.__H_0 = H_0                # [кДж/кг]   Располагаемый теплоперепад ступени (не полный)
                                  
        self.__fH_0         = None      # [кДж/кг]   Полный располагаемый теплоперепад ступени
        self.__fH_0s        = None      # [кДж/кг]   Полный располагаемый теплоперепад сопловой решетки
        self.__fH_0r        = None      # [кДж/кг]   Полный располагаемый теплоперепад Рабочей решетки
        
        # Углы ------------------------ # -----------------------------------------------------------------------------------
        
        self.__alpha_0 = alpha_0        # [deg]      Угол входа потока в ступень в абсолютном движении
        self.__alpha_1eef = alpha_1eef  # [deg]      Эффективный угол выхода потока из сопловую решетку в абсолютном движении 
        
        self.__alpha_1      = None      # [deg]      Угол выхода потока из сопловой решетки в абсолютном движении 
        self.__alpha_2      = None      # [deg]      Угол выхода потока из рабочей решетки в абсолютном движении
        self.__beta_1       = None      # [deg]      Угол входа потока в рабочую решетку в относительном движении
        self.__beta_2eef    = None      # [deg]      Эффективный угол выхода потока из рабочей решетки
        self.__beta_2       = None      # [deg]      Угол выхода потока из сопловой решетки в относительном движении

        self.__alpha_inst   = None      # [deg]      Установочный угол сопловой решетки
        self.__beta_inst    = None      # [deg]      Установочный угол рабочей решетки 
        
        # Скорости -------------------- # -----------------------------------------------------------------------------------
        
        self.__u            = None      # [м/c]      Окружная скорость 
        self.__c_f          = None      # [м/c]      Фиктивная скорость
        
        self.__c_0 = c_0                # [м/c]      Абсолютная скорость потока на входе в ступень (в сопловую решетку)
                                 
        self.__c_1t         = None      # [м/с]      Абсолютная теоретическая скорость на выходе из сопловой решетки
        self.__c_1          = None      # [м/c]      Абсолютная скорость потока на выходе из сопловой решетки
                                                                            
        self.__w_1t         = None      # [м/c]      Относительная теоретическая скорость на входе в сопловую решетку
        self.__w_1          = None      # [м/c]      Относительная скорость потока на выходе из сопловой решетки
                                        
        self.__w_2t         = None      # [м/с]      Относительная теоретическая скорость рабоей решетки
        self.__w_2          = None      # [м/c]      Относительная скорость потока на выходе из рабочей решетки 
                                        
        self.__c_2          = None      # [м/c]      Абсолютная скорость потока на выходе из рабочей решетки
                                        
        # Числа Маха ------------------ # ----------------------------------------------------------------------------
        
        self.__M_1t         = None      # []         Число Маха на выходе из сопловой решетке 
        self.__M_2t         = None      # []         Число Маха на выходе из рабочей решетки 
        
        # Геометрия лопаток ----------- # ----------------------------------------------------------------------------
        
        self.__Delta_pr = Delta_pr      # [м]        Перекрыша
        
        self.__l_1          = None      # [м]        Высота сопловой лопатки
        self.__l_2          = None      # [м]        Высота рабочей лопатки     
        
        self.__b_1          = None      # [м]        Хорда сопловой решетки
        self.__b_2          = None      # [м]        Хорда рабочей решетки
        
        self.__partiality   = None      # []         Степень парциальности
        
        self.__F_1          = None      # [м^2]      Площадь сопловой решетки
        self.__F_2          = None      # [м^2]      Площадь рабочей решетки
        
        self.__z_1          = None      # []         Число сопловых лопаток
        self.__z_2          = None      # []         Число рабочих лопаток
        
        self.__t_1          = None      # [м]        Шаг сопловой решетки
        self.__t_2          = None      # [м]        Шаг рабочей решетки
        
        self.__t_1_rel      = None      # []         Относительный шаг сопловой решетки
        self.__t_2_rel      = None      # []         Относительный шаг рабочей решетки
        
        self.__scale_1      = None      # []         Масштаб сопловой решетки, относительно модельной
        self.__scale_2      = None      # []         Масштаб рабочей решетки, относительно модельной
        
        self.__B_1          = None      # [м]        Ширина сопловой решетки
        self.__B_2          = None      # [м]        Ширина сопловой решетки
        
        self.__throat_1     = None      # [м]        Длина горла сопловой решетки
        self.__throat_2     = None      # [м]        Длина горла рабочей решетки
        
        self.__R_edge_in_1  = None      # [м]        Радиус входной кромки сопловой решетки
        self.__R_edge_out_1 = None      # [м]        Радиус выходной кромки сопловой решетки
        
        self.__R_edge_in_2  = None      # [м]        Радиус входной кромки сопловой решетки
        self.__R_edge_out_2 = None      # [м]        Радиус выходной кромки сопловой решетки
        
        
        
        
        # Потери ---------------------- # ----------------------------------------------------------------------------
        
        self.__Delta_Hs     = None      # [кДж/кг]   Потери в сопловой решетке
        self.__Delta_Hr     = None      # [кДж/кг]   Потери в рабочей решетке
        self.__Delta_Hout   = None      # [кДж/кг]   Потери с выходной скоростью
        self.__xi_sprof     = None      # []         Относительные профильные потери в сопловой решетке
        self.__xi_rprof     = None      # []         Относительные профильные потери в рабочей решетке
        self.__xi_ssec      = None      # []         Относительные суммарный потери в сопловой решетке
        self.__xi_rsec      = None      # []         Относительные суммарный потери в рабочей решетке
        self.__xi_sann      = None      # []         Относительные концевые потери в сопловой решетке
        self.__xi_rann      = None      # []         Относительные концевые потери в рабочей решетке 
        
        self.__eff_ol       = None      # []         Относительный лопаточный КПД ступени
        
        # Экономически-динамические п.  # ----------------------------------------------------------------------------
        
        self.__Ru           = None      # [кH]       Радиальное усилие
        self.__Lu           = None      # [кДж/кг]   Удельная работа ступени
        self.__Nu           = None      # [МВт]      Мощность ступени
        self.__E_0          = None      # [кДж/кг]   Располагаемая энергия ступени
        
        # Прочие величины ------------- # ----------------------------------------------------------------------------
        
        self.__G_0 = G_0                # [кг/s]     Массовый расход в ступени
        self.__d_hub = d_hub            # [м]        Корневой диаметр 
        self.__n = n                    # [1/s]      Частота вращения вала
        self.__reaction = reaction      # []         Степень рективности ступени на среднем диаметре
        self.__kappa_vs = kappa_vs      # []         Коэффициент выходной скорости
        
        self.__vane_grid    = None      # ...        Параметры сопловой решетки
        self.__rotor_grid   = None      # ...        Параметры рабочей решетки
        
        self.__isOK         = None      # ...        Просчитывается ли ступень
        
        # Основной расчет ступени 
        self.calculate()
        
    
    
    def __calculate_full_heat_transfer(self):
        # Полный располагаемый теплоперепад на ступень
        self.__fH_0 = self.__H_0 + np.power(self.__c_0, 2)/2000
        
        # Фиктивная скорость ступени
        self.__c_f = np.sqrt(2000*self.__fH_0)
        
    def __calculate_vane_thermal_points(self, phi):
        """
            Раасчитывает термодинамические точки на h,s-диаграмме для сопловой решетки
        """
        
        # Т.Т полных параметров на входе в ступень
        fh_0 = self.__Point_0.h(0) + np.power(self.__c_0, 2)/2000
        fs_0 = self.__Point_0.s(0)
        self.__fPoint_0 = ThPoint(h=fh_0,s=fs_0)
        
        # Располагаемый теплоперепад на сопловую решетку
        self.__fH_0s = self.__fH_0 * (1 - self.__reaction)
        
        # Т.Т статических теоретических параметров на выходе из сопловой решетки
        h_1t = self.__fPoint_0.h(0) - self.__fH_0s
        s_1t = self.__fPoint_0.s(0)
        self.__Point_1t = ThPoint(h=h_1t,s=s_1t)

        # Теоретическая скорость в сопловой решетке
        self.__c_1t = np.sqrt(2000*self.__fH_0s)
        self.__c_1 = self.__c_1t * phi
        
        # Расчет потерь в сопловой решетке
        self.__Delta_Hs = np.power(self.__c_1t, 2) / 2000 * (1 - np.power(phi, 2))
        
        # Т.Т. действительная статических параметром на входе в рабочую решетку
        h_1 = self.__Point_1t.h(0) + self.__Delta_Hs
        p_1 = self.__Point_1t.p(0)
        self.__Point_1 = ThPoint(p=p_1, h=h_1)
        
    def __calculate_rotor_thermal_points(self, psi):
        """
            Рассчитывает термодинамические точки на h,s-диаграмме для рабочей решетки
        """
        # Т.Т. полных параметров на входе в рабочую решетку
        fh_1 = self.__Point_1.h(0) + np.power(self.__w_1, 2) / 2000
        fs_1 = self.__Point_1.s(0)
        self.__fPoint_1 = ThPoint(h=fh_1, s=fs_1)
        
        # Т.Т. (2't)
        h_2tt = self.__fPoint_0.h(0) - self.__H_0 - np.power(self.__c_0, 2) / 2000
        s_2tt = self.__fPoint_0.s(0)
        self.__Point_2tt = ThPoint(h=h_2tt, s=s_2tt)
        
        # Т.Т. теоретических статических параметров на выходе из рабочей решетки
        p_2 = self.__Point_2tt.p(0)
        s_2t = self.__Point_1.s(0)
        self.__Point_2t = ThPoint(p=p_2, s=s_2t)
        
        # Полный теплоперепад на рабочую решетку
        self.__fH_0r = self.__fPoint_1.h(0) - self.__Point_2t.h(0)
        
        # Относительная теоретическая скорость в рабочей решетке
        self.__w_2t = np.sqrt(2000 * (self.__fH_0r))
        
        # Потери в рабочей решетке
        self.__Delta_Hr = np.power(self.__w_2t, 2) / 2000 * (1 - np.power(psi, 2))
        
        # Т.Т. действительных статических параметров на выходе из рабочей решетки
        h_2 = self.__Point_2t.h(0) + self.__Delta_Hr
        self.__Point_2 = ThPoint(p=p_2, h=h_2)   
            
    def __calculate_out_stage_thermal_points(self):
        """
            Рассчитывает термодинамические точки выхода из ступени
        """
        
        self.__X = self.__u / self.__c_f
        
        # Потери с выходной скоростью
        self.__Delta_Hout = np.power(self.__c_2, 2)/2000
        
        # Т.Т статических параметров на выходе из ступени
        p_outt = self.__Point_2.p(0)
        h_outt = self.__Point_2.h(0) + (1 - self.__kappa_vs) * self.__Delta_Hout
        self.__Point_outt = ThPoint(p=p_outt, h=h_outt)
        
        # Т.Т. полных парамеров на выходе из ступени
        s_out = self.__Point_outt.s(0)
        h_out = self.__Point_2.h(0) + self.__kappa_vs * self.__Delta_Hout
        self.__Point_out = ThPoint(h=h_out, s=s_out)
    
    def __find_critical_point(self, StartPoint:ThPoint, EndPoint:ThPoint) -> ThPoint:
        '''
            Для определения высоты лопатки при сверхкритическом течении необходимо
            определить критический объем и критическую скорость. Для этого 
            необходимо найти термодинамическую точку на h,s-диаграмме, соответствующую
            М = 1, а затем и остальные параметры в этой точке.
                    
            Для этого необходимо менять давления p_i =[fp_0...p_1], определять местную
            скорость звука a_i,h_i = f(p_i, fs_0) и c_it = sqrt(2*(fh_0 - h_i)), 
            далее просчитывать число Маха М, таким образом найдется число Маха, равное 
            единице, в этом сечении и будет критическая точка. 
                    
            Погрешность вычисления числа Маха возьмем Delta = 0.001
                    
            Так как число Маха в целом монотонно возрастает при убывании давления, то можно
            воспользоваться бинарным поиском. Однако непосредственно при приближении к Маху
            М = 1.00 +- 0.06 наблюдается не линейная зона изменения величины, поэтому
            целесообразно отказаться от бинарного поиска в пользу последовательного перебора.
            
            Для того, чтобы отследить зону, в которой число Маха перестает меняться (при при-
            ближении к 1), рассчитывается абслютная разность числа Маха на этой итерации и на 
            прошлой. Если число Маха перестает изменяться больше чем на 0.1, то необходимо на-
            чать перебор.
            
            Перебор осущесвляется с помощью умножения давления на 0.95 или 1.05 в зависимости
            от величины числа Маха (если М < 1, то умножаем на 0.95, если М > 1, то умножаем
            на 1.05).
            
            Важно, что при изначальном поиске давления (бинарном поиске) значение каждую итер-
            ацию изменяется на +- 0.1 от числа Маха, поэтому необходимо проверять разницу числа
            Маха для значений меньше 1 и больше - отдельно, так как иначе цикл никогда не за-
            вершится.
            
            Данное решение ПОКА ЧТО является самым эффективным.
        '''
        # Если точки не имеют одну энтропию
        if (np.abs(StartPoint.s(0) - EndPoint.s(0)) > 0.001):
            exc_values = "StartPoint.s={}, EndPoint.s={}".format(StartPoint.s(0), EndPoint.s(0))
            raise Exception(
                "Invalid input thermal points, entropy is not the same:{}".format(exc_values)
            )
        
        p_start = StartPoint.p(0)
        p_end = EndPoint.p(0)
        
        # Начальные приближения:
        p_i = np.mean([p_start, p_end])
        
        previous_M_higher = -1
        previous_M_lower  = -1
        # Поиск
        while (True):
            # i-ая термодинамическая точка 
            Point_i = ThPoint(p=p_i, s=StartPoint.s(0))
            # i-ая скорость звука
            a_i = Point_i.a(0)
            # i-ая энтальпия
            h_i = Point_i.h(0)
            # i-ая теоретическая абсолютная скорость
            c_it = np.sqrt(2*(StartPoint.h(0) - h_i)*1000)
            # i-ое число Маха
            M_i = c_it/a_i
            
            # Проверяем значение числа Маха
            # Если число Маха М ~ 1
            if (np.abs(M_i - 1) < 0.001):
                break
            
            # Среднее выполнение конструктора (не учитывая функции после данной) с разными
            # значениеми литералов:
            
            # 0.98,1.02 - ave time = 0.11454
            # 0.95,1.05 - ave time = 0.11394
            # 0.9, 1.1  - ave time = 0.11427
            
            # Если число Маха М > 1, то необходимо увеличить промежуточное давление
            elif M_i < 1:
                # если число Маха изменятся линейно
                if not(np.abs(previous_M_lower - M_i) < 0.1):
                    p_i = np.mean([p_i, p_end])
                # если число Маха перестало линейно изменяться
                else:
                    p_i = 0.98 * p_i
                previous_M_lower = M_i
                
            # Если число Маха М < 1, то необходимо уменьшить промежуточное давление
            elif M_i > 1:                
                # если число Маха изменятся линейно
                if not(np.abs(previous_M_higher - M_i) < 0.1):
                    p_i = np.mean([StartPoint.p(0), p_i])
                # если число Маха перестало линейно изменяться
                else:
                    p_i = 1.02 * p_i
                previous_M_higher = M_i
        
        return Point_i
        
    def __find_critical_point_1(self) -> ThPoint:
        # Начальные приближения:
        p_i = np.mean([self.__fp_0, self.__p_1])
        
        # Поиск
        while (True):
            # i-ая термодинамическая точка 
            Point_i = ThPoint(p=p_i, s=self.__fs_0)
            # i-ая скорость звука
            a_i = Point_i.a(0)
            # i-ая энтальпия
            h_i = Point_i.h(0)
            # i-ая теоретическая абсолютная скорость
            c_it = np.sqrt(2*(self.__fh_0 - h_i)*1000)
            # i-ое число Маха
            M_i = c_it/a_i
            
            # Проверяем значение числа Маха
            # Если число Маха М ~ 1
            if (np.abs(M_i - 1) < 0.001):
                break
            # Если число Маха М > 1, то необходимо увеличить промежуточное давление
            elif M_i < 1:
                if (np.abs(1 - M_i) > 0.100):
                    p_i = np.mean([p_i, self.__p_1])
                else:
                    p_i = 0.99 * p_i
                
            # Если число Маха М < 1, то необходимо уменьшить промежуточное давление
            elif M_i > 1:
                if (np.abs(1 - M_i) > 0.100):
                    p_i = np.mean([self.__fp_0, p_i])
                else:
                    p_i = 1.01 * p_i
        
        return Point_i
         
    def __get_normal_vane_area_and_alpha_1(self,phi,mu_1):
        """
            Функция предназначена для расчета площади сопловой решетки и угла выхода при числе Маха М < 1.
        """ 
        # Абсолютный угол выхода потока из сопловой решетки, в градусах
        # alpha_1 = np.rad2deg(np.arcsin(mu_1 * np.sin(np.deg2rad(self.__alpha_1eef))/phi))
        alpha_1 = self.__alpha_1eef
        
        # Выходная площадь сопловой решетки
        F_1 = self.__G_0 * self.__Point_1t.v(0) / mu_1 / self.__c_1t
    
        return F_1, alpha_1
    
    def __get_normal_rotor_area_and_rotor_angles(self, psi, mu_2):
        # Площадь рабочей решетки
        F_2 = self.__G_0 * self.__Point_2t.v(0) / mu_2 / self.__w_2t
            
        # Эффективный относительный угол выхода потока из рабочей решетки
        beta_2eef = np.rad2deg(np.arcsin(F_2 / (np.pi * self.__d * self.__l_2)))
            
        # Реальный относительный угол выхода потока из рабочей решетки
        beta_2 = np.rad2deg(np.arcsin(mu_2 * np.sin(np.deg2rad(beta_2eef)) / psi))
        
        return F_2, beta_2, beta_2eef
    
    def __get_critical_rotor_area_and_rotor_angles(self, psi, mu_2):
        # Поиск критической точки в сопловой решетке
        Point_cr = self.__find_critical_point(self.__fPoint_1, self.__Point_2t)

        # Критический удельный объем
        v_2cr = Point_cr.v(0)
                
        # Критическая абсолютная скорость в сопловой решетке
        w_2cr = np.sqrt(2000*(self.__fPoint_1.h(0)- Point_cr.h(0))) 
        
        # Приведенная относительная площадь
        q = self.__w_2t * v_2cr / w_2cr / self.__Point_1t.v(0)
        
        # Выходная площадь сопловой решетки
        F_2 = self.__G_0 * v_2cr / mu_2 / w_2cr
              
        # Эффективный относительный угол выхода потока из рабоей решетки, в градусах
        beta_2eef = np.rad2deg(np.arcsin(self.__F_2 / (np.pi * self.__d * self.__l_2)))
        
        # Реальный относительный угол выхода потока из рабочей решетки
        beta_2 = np.rad2deg(np.arcsin(mu_2 * np.sin(self.__beta_2eef) / psi / q))
        
        return F_2, beta_2, beta_2eef
        
    def __get_critical_vane_area_and_alpha_1(self,phi,mu_1):
        """
            Функция предназначена для расчета площади сопловой решетки и угла выхода при числе Маха М > 1.
        """ 
        # Поиск критической точки в сопловой решетке
        Point_cr = self.__find_critical_point(self.__fPoint_0, self.__Point_1t)

        # Критический удельный объем
        v_1cr = Point_cr.v(0)
                
        # Критическая абсолютная скорость в сопловой решетке
        c_1cr = np.sqrt(2000*(self.__fPoint_0.h(0)- Point_cr.h(0))) 
        
        # Приведенная относительная площадь
        q = self.__c_1t * v_1cr / c_1cr / self.__Point_1t.v(0)
                
        # Абсолютный угол выхода потока из сопловой решетки, в градусах
        # alpha_1 = np.rad2deg((np.arcsin(mu_1 * np.sin(np.deg2rad(self.__alpha_1eef)/phi/q))))
        alpha_1 = np.rad2deg((np.arcsin(np.sin(np.deg2rad(self.__alpha_1eef)/q))))
        
        # Выходная площадь сопловой решетки
        F_1 = self.__G_0 * v_1cr / mu_1 / c_1cr
        
        return F_1, alpha_1
        
    def __calculate_vane_area_and_length(self, phi, mu_1):
        '''
            Функция рассчитывает длину сопловых лопаток на основе коэффициентов скорости и расхода.
              
        '''
        # Число маха для сопловой решетки
        self.__M_1t = self.__c_1t/self.__Point_1t.a(0)
    
        # Выходная площадь сопловой решетки
        # Если число Маха M < 1
        if (self.__M_1t < 1):
            # Выходная площадь сопловой решетки & Абсолютный угол выхода потока из сопловой решетки
            self.__F_1, self.__alpha_1 = self.__get_normal_vane_area_and_alpha_1(phi=phi, mu_1 = mu_1)
                
        # Если число Маха М > 1
        else:
            # Выходная площадь сопловой решетки & Абсолютный угол выхода потока из сопловой решетки
            self.__F_1, self.__alpha_1 = self.__get_critical_vane_area_and_alpha_1(phi=phi, mu_1=mu_1)
        
        # Коэффициенты квадратного уравнения:
        koef_a = np.pi * np.sin(np.deg2rad(self.__alpha_1))
        koef_b = np.pi * self.__d_hub * np.sin(np.deg2rad(self.__alpha_1))
        
        # Выходная высота сопловых лопаток 
        self.__l_1 = (-koef_b + np.sqrt(koef_b**2 + 4*koef_a*self.__F_1))/(2 * koef_a)
        
        # Средний диаметр 
        self.__d = self.__d_hub + self.__l_1
        
        # Окружная скорость
        self.__u = np.pi * self.__d * self.__n
        
        # Если лопатки получились меньше чем 0.013 мм
        if (self.__l_1 < 0.013):
            # Степень парциальности != 1
            self.__partiality = 6 * np.sqrt(self.__l_1)
            self.__l_1 = self.__l_1 / self.__partiality
        else:
            # Степень парциальности = 1
            self.__partiality = 1
    
    def __calculate_vane_speed_triangle(self, phi):
        """
            Функция рассчитывает величины соплового треугольника скоростей.
            Список рассчитываемых величины:
            с_1         - Действительная абсолютная скорость на выходе из сопловой решетки
            beta_1      - Угол выхода из сопловой решетки в относительном движении
            w_1t,w_1    - Относительная теоретическая и действительная скорость на входе в рабочую решетку
        """
        # Действительная абсолютная скорость на выходе из сопловой решетки
        self.__c_1 = self.__c_1t * phi
        
        # Действительная относительная скорость на входе в рабочую решетку
        self.__w_1 = np.sqrt(np.power(self.__c_1, 2) + np.power(self.__u, 2) - 2 * self.__u * self.__c_1 * np.cos(np.deg2rad(self.__alpha_1)))
        
        # Теоретическая относительная скорость на входе в рабочую решетку
        # self.__w_1t = self.__w_1 / psi
        
        # Угол выхода из сопловой решетки в относительном движении
        self.__beta_1 = self.__alpha_1 + np.rad2deg(np.arccos((np.power(self.__c_1, 2) + np.power(self.__w_1, 2) - np.power(self.__u, 2))/(2*self.__c_1*self.__w_1))) 
    
    def __calculate_rotor_area_and_length(self, psi, mu_2):
        """
            Рассчет длины и площади рабочей лопатки
        """
        # Длина рабочей лопатки (с учетом перекрыши)
        self.__l_2 = self.__l_1 + self.__Delta_pr
        
        # Число маха на выходе из рабочей решетки
        self.__M_2t = self.__w_2t / self.__Point_2t.a(0)
        
        # Далее расчет изменяется для решеток с сверхзвуком и дозвуком
        # Если дозвуковое течение
        if (self.__M_2t < 1):
            # Площадь рабочей лопатки, реальный и эффективный угол выхода потока
            self.__F_2, self.__beta_2, self.__beta_2eef = self.__get_normal_rotor_area_and_rotor_angles(psi=psi,mu_2=mu_2)   
        # Если сверхзвуковое течение
        else:
            # Площадь рабочей лопатки, реальный и эффективный угол выхода потока
            self.__F_2, self.__beta_2, self.__beta_2eef = self.__get_critical_rotor_area_and_rotor_angles(psi=psi,mu_2=mu_2)
    
    def __calculate_rotor_speed_triangle(self, psi):
        # Реальная относительная скорость в рабочей решетке
        self.__w_2 = self.__w_2t * psi
        
        # Действительная абсолютная скорость на выходе из рабочей решетки
        self.__c_2 = np.sqrt(np.power(self.__w_2, 2) + np.power(self.__u, 2) - 2 * self.__u * self.__w_2 * np.cos(np.deg2rad(self.__beta_2)))
        
        # Абсолютный выход потока из рабоче решетки
        self.__alpha_2 = self.__beta_2 + np.rad2deg(np.arccos((np.power(self.__c_2,2)+np.power(self.__w_2,2)-np.power(self.__u,2))/(2*self.__w_2*self.__c_2)))
    
    def __select_vane_grid(self, b_1):
        """
            Выбирает сопловую решетку из таблицы, расположенной в Щегляеве, том 1
        """
        self.__vane_grid = GridProfile(
            M=self.__M_1t,           # Число Маха в рабочей решетке
            type=0,                  # Тип решетки, 1 - рабочая
            in_angle=self.__alpha_0, # Входной угол
            out_angle=self.__alpha_1 # Выходной угол
        )
        # Если решетка не подобралась -> исключение
        if not self.__vane_grid.isOK():
            raise Exception('Vane Grid was not selected 0_0')
        
        # Дробное число сопловых лопаток
        self.__z_1 = np.pi * self.__d / self.__vane_grid.get_rel_t() / b_1
        
        # Целое число сопловых лопаток (четное)
        self.__z_1 = round(self.__z_1 / 2) * 2
        
        # Шаг сопловой решетки
        self.__t_1 = np.pi * self.__d / self.__z_1
        
        # Относительный шаг сопловой решетки
        self.__t_1_rel = self.__t_1 / b_1
        
        # Рассчет установочного угла
        self.__alpha_inst = self.__vane_grid.calculate_inst_angle(in_angle=self.__alpha_0, t_rel = self.__t_1_rel)
        
        # Масштаб для сопловой решетки
        self.__scale_1 = b_1 / self.__vane_grid.get_b()
        
        # Ширина сопловой решетки
        self.__B_1 = b_1 * np.sin(np.deg2rad(self.__alpha_inst))
        
        # Горло сопловой решетки
        self.__throat_1 = self.__scale_1 * self.__vane_grid.get_a()
        
        # Радиусы входной и выходной кромок
        self.__R_edge_in_1 = self.__scale_1 * self.__vane_grid.get_Delta_in()
        self.__R_edge_out_1 = self.__scale_1 * self.__vane_grid.get_Delta_out()
        
        
    
    def __select_rotor_grid(self,b_2):
        """
            Выбирает рабочую решетку из таблицы, расположенной в Щегляеве, том 1
        """
        
        self.__rotor_grid = GridProfile(
            M=self.__M_2t,          # Число Маха в рабочей решетке
            type=1,                 # Тип решетки, 1 - рабочая
            in_angle=self.__beta_1, # Входной угол
            out_angle=self.__beta_2 # Выходной угол
        )
        # Если решетка не подобралась -> исключение
        if not self.__rotor_grid.isOK():
            raise Exception('Rotor Grid was not selected 0_0')
        
        # Дробное число рабочих лопаток
        self.__z_2 = np.pi * self.__d / self.__rotor_grid.get_rel_t() / b_2
        
        # Целое число рабочих лопаток (четное)
        self.__z_2 = int(self.__z_1)
        
        # Шаг рабочих решетки
        self.__t_2 = np.pi * self.__d / self.__z_2
        
        # Относительный шаг рабочих решетки
        self.__t_2_rel = self.__t_2 / b_2
        
        # Рассчет установочного угла
        self.__beta_inst = self.__rotor_grid.calculate_inst_angle(in_angle=self.__beta_1, t_rel = self.__t_2_rel)
        
        # Масштаб для сопловой решетки
        self.__scale_2 = b_2 / self.__rotor_grid.get_b()
        
        # Ширина сопловой решетки
        self.__B_2 = b_2 * np.sin(np.deg2rad(self.__beta_inst))
        
        # Горло сопловой решетки
        self.__throat_2 = self.__scale_2 * self.__rotor_grid.get_a()
        
        # Радиусы входной и выходной кромок
        self.__R_edge_in_2 = self.__scale_2 * self.__rotor_grid.get_Delta_in()
        self.__R_edge_out_2 = self.__scale_2 * self.__rotor_grid.get_Delta_out()
        
        
    # TODO
    def __check_bending_stresses(self, b_2) -> bool:
        # Момент сопротивления по кромкам по оси XX
        W_min = self.__rotor_grid.get_Wxx_edge() * np.power(b_2 / self.__rotor_grid.get_b(), 3)
        
        # Число рабочих лопаток
        z_2 = np.pi * self.__d / b_2 / self.__rotor_grid.get_rel_t()
        
        # Целое число рабочих лопаток
        if (abs(int(z_2) - z_2) != 0):
            z_2new = int(z_2) + 1
        else:
            z_2new = int(z_2)
            
        # Радиальное усилие
        Ru = self.__G_0 * (self.__w_1 * np.cos(np.deg2rad(self.__beta_1)) + self.__w_2 * np.cos(np.deg2rad(self.__beta_2)))
        
        # Изгибающее напряжение
        sigma_bend = Ru * self.__l_2 / 2 / z_2new / W_min / 1_000_000
        
        # TODO проверка только для стали
        if (sigma_bend < 27.5):
            return True
        else:
            return False

    def __calculate_losses(self, type, b):
        '''
        Расчет потерь
        '''
    
        if (type == 's'):
            self.__xi_sprof, self.__xi_ssec, self.__xi_sann = calculate_losses(
                inlet_angle = self.__alpha_0,
                outlet_angle = self.__alpha_1,
                blade_inlet_angle = self.__alpha_0 - 5,
                design_inc = 0,
                axial_chord = self.__B_1,
                throat = self.__throat_1,
                pitch = self.__t_1,
                te_radius = self.__R_edge_out_1,
                chord = b,
                height = self.__l_1,
                nu = (self.__Point_0.kv() + self.__Point_1t.kv())/2,
                ks = 1e-5,
                w1 = self.__c_0,
                w2 = self.__c_1
            )
            self.__xi_sann = self.__xi_ssec - self.__xi_sprof
            
            return self.__xi_ssec
        elif (type == 'r'):
            self.__xi_rprof, self.__xi_rsec, self.__xi_rann = calculate_losses(
                inlet_angle = self.__beta_1,
                outlet_angle = self.__beta_2,
                blade_inlet_angle = self.__beta_1 - 5,
                design_inc = 0,
                axial_chord = self.__B_2,
                throat = self.__throat_2,
                pitch = self.__t_2,
                te_radius = self.__R_edge_out_2,
                chord = b,
                height = self.__l_2,
                nu = (self.__Point_1.kv() + self.__Point_2t.kv())/2,
                ks = 1e-5,
                w1 = self.__w_1,
                w2 = self.__w_2
            )
            self.__xi_rann = self.__xi_rsec - self.__xi_rprof
            
            return self.__xi_rsec
        else:
            raise Exception('Unknown grid type in input parameters')
        
    
    def __calculate_power_efficiency_radial_force(self):
        """Функция рассчитывает мощность, КПД, радиальные усилия в ступени
        """
        # Окружное усилие определим по нескольким формулам для проверки
        alpha_1_rad = np.deg2rad(self.__alpha_1)
        alpha_2_rad = np.deg2rad(self.__alpha_2)
        beta_1_rad  = np.deg2rad(self.__beta_1)
        beta_2_rad  = np.deg2rad(self.__beta_2)
        
        Ru_1 = self.__G_0 * (self.__c_1 * np.cos(alpha_1_rad) + self.__c_2 * np.cos(alpha_2_rad)) / 1000
        Ru_2 = self.__G_0 * (self.__w_1 * np.cos(beta_1_rad) + self.__w_2 * np.cos(beta_2_rad)) / 1000
        
        # Если не совпали, то выкидываем исключение
        if (round(Ru_1, 3) != round(Ru_2, 3)):
            raise Exception("Invalid radial force calculation, Ru_1 != Ru_2 -> {} != {}".format(Ru_1, Ru_2))
        else:
            self.__Ru = Ru_1
        
        # Удельная работа ступени
        Lu_1 = self.__Ru * self.__u / self.__G_0
        Lu_2 = (np.power(self.__c_1, 2) - np.power(self.__c_2, 2) + np.power(self.__w_2, 2) - np.power(self.__w_1, 2))/2000
        Lu_3 = self.__fH_0 - self.__Delta_Hr - self.__Delta_Hs - self.__Delta_Hout
        
        if ((round(Lu_1, 3) != round(Lu_2, 3)) &
            (round(Lu_2, 3) != round(Lu_3, 3)) &
            (round(Lu_3, 3) != round(Lu_1, 3))):
            raise Exception("Invalid specific stage work calculation, Lu_1 != Lu_2 != Lu_1 -> {} != {} != {}".format(Lu_1, Lu_2, Lu_3))
        else:
            self.__Lu = Lu_3
        
        # Мощность ступени
        self.__Nu = self.__Lu * self.__G_0 / 1000
        
        # Располагаемая энергия ступени
        self.__E_0 = self.__fH_0 - self.__kappa_vs * self.__Delta_Hout
        
        # Относительный лопаточный КПД
        eff_ol_1 = self.__Lu / self.__E_0
        eff_ol_2 = (self.__E_0 - self.__Delta_Hs - self.__Delta_Hr - (1 - self.__kappa_vs) * self.__Delta_Hout)/self.__E_0
        
        if (round(eff_ol_1, 3) != round(eff_ol_2, 3)):
            raise Exception("Invalid efficency calculation, eff_ol_1 != eff_ol_2 -> {} != {}".format(eff_ol_1, eff_ol_2))   
        else:
            self.__eff_ol = eff_ol_1
    
    def calculate(self):
        '''
        Последовательный расчет ступени
        '''
        
        # 1. Первое минимальное приближение хорд лопаток
        b_1 = 0.030
        b_2 = 0.020
       
        phi = 0.969
        mu_1 = 0.975
        
        psi = 0.952
        mu_2 = 0.959
        
        isBendingOK = False
        isSpeedKoefsOK = False
        
        # 2. Рассчитываем полный теплоперепад в ступени
        self.__calculate_full_heat_transfer()
        
        # 3.
        while not(isBendingOK * isSpeedKoefsOK):
            
            # 3.1. Рассчитываем параметры для сопловой решетки в h,s-диаграммы
            self.__calculate_vane_thermal_points(phi=phi)
            
            # 3.2. Рассчитываем площадь и длину сопловой лопатки
            self.__calculate_vane_area_and_length(phi=phi, mu_1=mu_1)
            
            # 3.3. Рассчитываем скорости и углы для сопловой решетки из треугольника скоростей
            self.__calculate_vane_speed_triangle(phi=phi)
            
            
            # 3.4. Рассчитываем параметры для рабочей решетки в h,s-диаграмме
            self.__calculate_rotor_thermal_points(psi=psi)
            
            # 3.5. Рассчитываем площадь и длину рабочей лопатки
            self.__calculate_rotor_area_and_length(psi=psi, mu_2=mu_2)
            
            # 3.6.0 Подбираем сопловую решетку
            self.__select_vane_grid(b_1=b_1)
            
            # 3.6.1 Подбираем рабочую решетку
            self.__select_rotor_grid(b_2=b_2)
            
            # 3.7. Рассчитываем треугольник скоростей в рабочей решетке 
            self.__calculate_rotor_speed_triangle(psi=psi)
            
            # 3.8. Проверяем рабочую решетку на изгиб
            if not(self.__check_bending_stresses(b_2=b_2)):
                b_1 += 0.005
                b_2 += 0.005
                continue
            else:
                isBendingOK = True
                
            # 3.9. Если прошлое условие на изгиб выполнилось, начинаем уточнять коэффициенты скорости
            # Для сопловой решетки 
            new_phi = np.sqrt(1 - self.__calculate_losses(type='s',b=b_1))
            # Для рабочей решетки
            new_psi = np.sqrt(1 - self.__calculate_losses(type='r',b=b_2))
            
            # Если совпали коэффициенты
            if ((round(new_phi, 3) == round(phi, 3)) & (round(new_psi, 3) == round(psi, 3))):
                isSpeedKoefsOK = True
                
            # Если коэффициенты не совпали, то совершаем расчет заново
            if (round(new_phi, 3) != round(phi, 3)):
                phi = new_phi
            if (round(new_psi, 3) != round(psi, 3)):
                psi = new_psi
            
            
        self.__b_1 = b_1
        self.__b_2 = b_2    
        self.__phi = phi
        self.__psi = psi
        
        
        # 4. Расчет т.т на выходе ступени
        self.__calculate_out_stage_thermal_points()
        
        # 5. Расчет КПД, мощности, уточнение усилий
        self.__calculate_power_efficiency_radial_force()
        

    def get_results(self):
        return {
            'd':        self.__d,
            'l_1':      self.__l_1, 
            'l_2':      self.__l_2, 
            'c_0':      self.__c_0, 
            'c_1':      self.__c_1, 
            'c_2':      self.__c_2, 
            'w_1':      self.__w_1,
            'w_2':      self.__w_2,
            'alpha_0':  self.__alpha_0,
            'alpha_1':  self.__alpha_1,
            'beta_1':   self.__beta_1,
            'beta_2':   self.__beta_2,
            'u':        self.__u,
            'X':        self.__X,
            'Delta_Hs': self.__Delta_Hs,
            'Delta_Hr': self.__Delta_Hr,
            'phi':      self.__phi,
            'psi':      self.__psi,
            'eff_ol':   self.__eff_ol,
            'E_0':      self.__E_0,
            'b_1':      self.__b_1,
            'b_2':      self.__b_2,
            'Name_1':   self.__vane_grid.get_name(),
            'Name_2':   self.__rotor_grid.get_name()
        }
    
    def get_l_1(self): return self.__l_1
    def get_reaction(self): return self.__reaction
    def get_d(self): return self.__d

class ExistingGrid:
    """
        Класс создан для хранения основных параметров определенной решетки.
    """
    
    def __init__(       
        self,
        # Геометрические характеристики
        # ----------------------------------- # ---------------------------------------------------------------
        b:float,                              # [м]      Хорда 
        B:float,                              # [м]      Ширина
        min_sec_length:float,                 # [м]      Горло
        rel_t:float,                          # []       Относительный шаг
        ins_angle:float,                      # [deg]    Установочный угол
        in_edge:float,                        # [м]      Входная кромка
        out_edge:float,                       # [м]      Выходная кромка
        inst_angle_expr:str,                  #          Выражение для определения установочного угла
        #------------------------------------ # ---------------------------------------------------------------
        # Прочностные характеристики   
        # ----------------------------------- # ---------------------------------------------------------------
        x_weight_center:float,                # [м]      x-координата центра тяжести
        y_weight_center:float,                # [м]      y-координата центра тяжестт
        profile_area:float,                   # [м^2]    Площадь профиля
        inertia_moment_xx:float,              # [м^4]    Момент инерции относительно ХХ
        inertia_moment_yy:float,              # [м^4]    Момент инерции относительно YY
        resistance_moment_xx_back:float,      # [м^3]    Момент сопротивления относительно ХХ, по спинке
        resistance_moment_xx_edge:float,      # [м^3]    Момент сопротивления относительно ХХ, по кромкам
        resistance_moment_yy_in_edge:float,   # [м^3]    Момент сопротивления относительно YY, по вх. кромке
        resistance_moment_yy_out_edge:float,  # [м^3]    Момент сопротивления относительно YY, по вых. кромке
        
        ):
        
        # КонструктоРРР
        self.__b = b
        self.__B = B
        self.__min_sec_length = min_sec_length
        self.__rel_t = rel_t
        self.__ins_angle = ins_angle
        self.__in_edge = in_edge
        self.__out_edge = out_edge
        
        self.__x_weight_center = x_weight_center
        self.__y_weight_center = y_weight_center
        self.__profile_area = profile_area
        self.__inertia_moment_xx = inertia_moment_xx
        self.__inertia_moment_yy = inertia_moment_yy
        self.__resistance_moment_xx_back = resistance_moment_xx_back
        self.__resistance_moment_xx_edge = resistance_moment_xx_edge
        self.__resistance_moment_yy_in_edge = resistance_moment_yy_in_edge
        self.__resistance_moment_yy_out_edge = resistance_moment_yy_out_edge
        self.__inst_angle_expr = inst_angle_expr
        
    def __str__(self) -> str:
        out_str = ""
        out_str += 'b = {}\n'.format(self.__b)
        out_str += 'B = {}\n'.format(self.__B)
        out_str += 'min_sec_length = {}\n'.format(self.__min_sec_length)
        out_str += 'rel_t = {}\n'.format(self.__rel_t)
        out_str += 'ins_angle = {}\n'.format(self.__ins_angle)
        out_str += 'in_edge = {}\n'.format(self.__in_edge)
        out_str += 'out_edge = {}\n'.format(self.__out_edge)
        out_str += 'x_weight_center = {}\n'.format(self.__x_weight_center)
        out_str += 'y_weight_center = {}\n'.format(self.__y_weight_center)
        out_str += 'profile_area = {}\n'.format(self.__profile_area)
        out_str += 'inertia_moment_xx = {}\n'.format(self.__inertia_moment_xx)
        out_str += 'inertia_moment_yy = {}\n'.format(self.__inertia_moment_yy)
        out_str += 'resistance_moment_xx_back = {}\n'.format(self.__resistance_moment_xx_back)
        out_str += 'resistance_moment_xx_edge = {}\n'.format(self.__resistance_moment_xx_edge)
        out_str += 'resistance_moment_yy_in_edge = {}\n'.format(self.__resistance_moment_yy_in_edge)
        out_str += 'resistance_moment_yy_out_edge = {}\n'.format(self.__resistance_moment_yy_out_edge)
        out_str += 'inst_angle_expr = {}\n'.format(self.__inst_angle_expr)
        return out_str
    
    # Getters      
    def get_b(self): return self.__b
    def get_B(self): return self.__B
    def get_a(self): return self.__min_sec_length
    def get_rel_t(self): return self.__rel_t
    def get_standart_ins_angle(self): return self.__ins_angle
    def get_Delta_in(self): return self.__in_edge
    def get_Delta_out(self): return self.__out_edge
    def get_center_X(self): return self.__x_weight_center
    def get_center_Y(self): return self.__y_weight_center
    def get_F(self): return self.__profile_area
    def get_Ixx(self): return self.__inertia_moment_xx
    def get_Iyy(self): return self.__inertia_moment_yy
    def get_Wxx_back(self): return self.__resistance_moment_xx_back
    def get_Wxx_edge(self): return self.__resistance_moment_xx_edge
    def get_Wyy_in_edge(self): return self.__resistance_moment_yy_in_edge
    def get_Wyy_out_edge(self): return self.__resistance_moment_yy_out_edge
    def get_inst_angle_expr(self): return self.__inst_angle_expr

class GridProfile:
    """
        Класс, храняющий основные параметры турбинных сопловых и рабочих решеток МЭИ.
    """
    __data = pd.read_excel('API/grids.xlsx')
    
    def __init__(self, M:float=0.5, type=0, in_angle:float=90, out_angle:float=14, name:str = ""):
        try:
            # Копируем необходимую строку параметров для выбранной решетки
            if (name != ""):
                loc_par = self.__data[self.__data['Name'] == name]
                self.__name = name
            else:
                # Ограничение на включенные решетки (ручное отключение в таблице)
                loc_par = self.__data[self.__data['isOn'] == 1]
                
                # Выбор типа решетки: сопловая/рабочая
                loc_par = loc_par[loc_par['type'] == type]
                
                # Ограничение по числу Маха
                loc_par = loc_par[(loc_par['M_b'] <= M) & (loc_par['M_e'] > M)]
                
                # Ограничение по входному углу
                loc_par = loc_par[(loc_par['out_angle_b'] <= out_angle) & (loc_par['out_angle_e'] > out_angle)]
                
                # Ограничение по выходному углу
                loc_par = loc_par[(loc_par['in_angle_b'] <= in_angle) & (loc_par['in_angle_e'] > in_angle)]
                
            # self.__grids_list = loc_par
            
            # Обновляем индексы строк
            # self.__grids_list.reset_index(drop=True, inplace=True)
            loc_par.reset_index(drop=True, inplace=True)
            
            if loc_par.empty: self.__name_list = []
            else: self.__name_list = loc_par['Name'].tolist()
            
            
            # Если строка найдена не одна или не найдена
            # if (len(loc_par) != 1):
            #     raise Exception("Invalid number of founded grids, number of grids: {}".format(len(loc_par)))
            

            self.__name = (loc_par.loc[0, 'Name'])
            self.__existing_grid = ExistingGrid(
                b                               = float(loc_par.loc[0, 'b']),
                B                               = float(loc_par.loc[0, 'B']),
                min_sec_length                  = float(loc_par.loc[0, 'a']),
                rel_t                           = float(loc_par.loc[0, 'rel_t']),
                ins_angle                       = float(loc_par.loc[0, 'ins_angle']),
                in_edge                         = float(loc_par.loc[0, 'Delta_in']),
                out_edge                        = float(loc_par.loc[0, 'Delta_out']),
                x_weight_center                 = float(loc_par.loc[0, 'center_X']),
                y_weight_center                 = float(loc_par.loc[0, 'center_Y']),
                profile_area                    = float(loc_par.loc[0, 'F']),
                inertia_moment_xx               = float(loc_par.loc[0, 'Ixx']),
                inertia_moment_yy               = float(loc_par.loc[0, 'Iyy']),
                resistance_moment_xx_back       = float(loc_par.loc[0, 'Wxx_back']),
                resistance_moment_xx_edge       = float(loc_par.loc[0, 'Wxx_edge']),
                resistance_moment_yy_in_edge    = float(loc_par.loc[0, 'Wyy_in_edge']),
                resistance_moment_yy_out_edge   = float(loc_par.loc[0, 'Wyy_out_edge']),
                inst_angle_expr                 = str(loc_par.loc[0, 'ins_angle_expression'])
            )
            
        except Exception as ex:
            print(*ex.args)
    
    def __str__(self) -> str:
        loc_str = ''
        # loc_str += 'name = {}\n'.format(str(self.__name.to_list()[0]))
        loc_str += str(self.__name) + '\n'
        loc_str += str(self.__existing_grid)
        return loc_str
        
            
    # Getters   
    def get_name(self): return self.__name
    def get_name_list(self): return self.__name_list
    def get_grids_list(self): return self.__grids_list 
    def get_b(self): return self.__existing_grid.get_b()
    def get_B(self): return self.__existing_grid.get_B()
    def get_a(self): return self.__existing_grid.get_a()
    def get_rel_t(self): return self.__existing_grid.get_rel_t()
    def get_standart_ins_angle(self): return self.__existing_grid.get_standart_ins_angle()
    def get_Delta_in(self): return self.__existing_grid.get_Delta_in()
    def get_Delta_out(self): return self.__existing_grid.get_Delta_out()
    def get_center_X(self): return self.__existing_grid.get_center_X()
    def get_center_Y(self): return self.__existing_grid.get_center_Y()
    def get_F(self): return self.__existing_grid.get_F()
    def get_Ixx(self): return self.__existing_grid.get_Ixx()
    def get_Iyy(self): return self.__existing_grid.get_Iyy()
    def get_Wxx_back(self): return self.__existing_grid.get_Wxx_back()
    def get_Wxx_edge(self): return self.__existing_grid.get_Wxx_edge()
    def get_Wyy_in_edge(self): return self.__existing_grid.get_Wyy_in_edge()
    def get_Wyy_out_edge(self): return self.__existing_grid.get_Wyy_out_edge()
    
    def calculate_inst_angle(self, in_angle, t_rel, isAcc = True):
        if isAcc:
            if 'angle' in self.__existing_grid.get_inst_angle_expr():
                return eval(self.__existing_grid.get_inst_angle_expr().format(angle = in_angle, t_rel = t_rel))
            else:
                return float(self.__existing_grid.get_inst_angle_expr())
        else:
            return self.get_standart_ins_angle()
        
    def isOK(self):
        if len(self.__name_list) != 0: return True
        else: return False