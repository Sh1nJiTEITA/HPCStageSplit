# Используемые библиотеки

# Общие
import numpy as np
import pandas as pd

# Собственные
from . import ThPoint





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
        
        
        
        isStageOptimal:bool = True, # Рассчитывается ли ступень на оптимальные параметры
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
        
        # Входные величины
        #
        self.__p_0 = p_0                # [МПа]      Статическое давление на входе в ступень
        self.__t_0 = t_0                # [degC]     Статическая температура на входи в ступень
        self.__G_0 = G_0
        self.__d_hub = d_hub            # [м]        Корневой диаметр 
        self.__n = n                    # [1/s]      Частота вращения вала
        self.__reaction = reaction      # []         Степень рективности ступени на среднем диаметре
        self.__alpha_0 = alpha_0        # [deg]      Угол входа потока в ступень в абсолютном движении
        self.__alpha_1eef = alpha_1eef  # [deg]      Эффективный угол выхода потока из сопловую решетку в абсолютном движении 
        self.__H_0 = H_0                # [кДж/кг]   Располагаемый теплоперепад ступени (не полный)
        self.__c_0 = c_0                # [м/c]      Абсолютная скорость потока на входе в ступень (в сопловую решетку)
        self.__kappa_vs = kappa_vs      # []         Коэффициент выходной скорости
        self.__Delta_pr = Delta_pr      # [м]        Перекрыша
        
        self.__is_stage_optimal:bool = True # Флаг оптимальности ступени

        # Дополнительные параметры
        if not len(extra_param) == 0:
            for param in extra_param:
                # Если задано u/c_f, то ступень не оптимальная
                if (param == "u/c_f"):
                    self.__is_stage_optimal = False
                    self.__X = extra_param.get(param)
        
        
        # Выходные величины 
        #
        self.__l_1      = None          # [м]        Высота сопловой лопатки
        self.__l_2      = None          # [м]        Высота рабочей лопатки
        self.__c_1      = None          # [м/c]      Абсолютная скорость потока на выходе из сопловой решетки
        self.__c_2      = None          # [м/c]      Абсолютная скорость потока на выходе из рабочей решетки
        self.__w_1      = None          # [м/c]      Относительная скорость потока на выходе из сопловой решетки
        self.__w_2      = None          # [м/c]      Относительная скорость потока на выходе из рабочей решетки 
        self.__beta_1   = None          # [deg]      Угол входа потока в рабочую решетку в относительном движении
        self.__beta_2   = None          # [deg]      Угол выхода потока из сопловой решетки в относительном движении
        self.__u        = None          # [м/c]      Окружная скорость 
        #self.__X        = None          # []         Значение (u/c_f)_отп
        self.__Delta_Hs = None          # [кДж/кг]   Потери в сопловой решетке
        self.__Delta_Hr = None          # [кДж/кг]   Потери в рабочей решетке
        self.__xi_s     = None          # []         Относительные потери в сопловой решетке
        self.__xi_r     = None          # []         Относительные потери в рабочей решетке
        self.__eff_ol   = None          # []         Относительный лопаточный КПД ступени
        self.__E_0      = None          # [кДж/кг]   Располагаемая энергия ступени
        
        # Промежуточные величины
        #
        self.__phi      = None          # []         Коэффициент скорости для соплвоой решетки
        self.__psi      = None          # []         Коэффициент скорости для рабочей решетки
        self.__fH_0     = None          # [кДж/кг]   Полный располагаемый теплоперепад ступени
        self.__fH_0s    = None          # [кДж/кг]   Полный располагаемый теплоперепад сопловой решетки
        self.__c_1t     = None          # [м/с]      Абсолютная теоретическая скорость на выходе из сопловой решетки
        self.__w_1t     = None          # [м/c]      Относительная теоретическая скорость на входе в сопловую решетку
        self.__w_2t     = None          # [м/с]      Относительная теоретическая скорость рабоей решетки
        self.__fH_0r    = None          # [кДж/кг]   Полный располагаемый теплоперепад Рабочей решетки
        self.__p_1      = None          # [МПа]      Давление на выходе из сопловой решетки
        self.__p_2      = None          # [МПа]      Статическое давление на выходе из ступени
        self.__M_1t     = None          # []         Число Маха на выходе из сопловой решетке (в горле)
        self.__mu_1     = None          # []         Коэффициент расхода для сопловой решетки
        self.__mu_2     = None          # []         Коэффициент расхода для рабочей решетки
        self.__fp_0     = None          # [МПа]      Давление полного торможения на входе в ступень
        self.__fh_0     = None          # [кДж/кг]   Энтальпия полного торможения на входе в ступень
        self.__fs_0     = None          # [кДж/кг/К] Энтропия полного торможения на входе в ступень
        self.__alpha_1  = None          # [deg]      Абсолютный угол выхода потока из сопловой решетки
        self.__partiality = None        # []         Степень парциальности
        self.__c_f      = None          # [м/c]      Фиктивная скорость
        self.__v_1t     = None          # [м^3/кг]   Теоретический удельный объем после сопловой решетки
        self.__F_1      = None          # [м^2]      Площадь сопловой решетки
        self.__h_1t     = None          # [кДж/кг]   Теоретическая статическая энтальпия на входе в сопловую решетку
        self.__h_1      = None          # [кДж/кг]   Действительная статическая энтальпия на входе в сопловую решетку
        self.__h_2t     = None          # [кДж/кг]   Теоретическая энтальпия на выходе из рабочей решетки
        self.M_i        = None
         
        #self.__calculate_vane_length(phi=0.973, mu=0.978)
        self.__calculate_vane_length(phi=np.mean([0.93,0.96]), mu=np.mean([0.95,0.97]))
    
    
    def __find_critical_point_2(self) -> ThPoint:
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
        
        # Начальные приближения:
        p_i = np.mean([self.__fp_0, self.__p_1])
        
        previous_M_higher = self.__M_1t
        previous_M_lower  = self.__M_1t
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
            
            # Среднее выполнение конструктора (не учитывая функции после данной) с разными
            # значениеми литералов:
            
            # 0.98,1.02 - ave time = 0.11454
            # 0.95,1.05 - ave time = 0.11394
            # 0.9, 1.1  - ave time = 0.11427
            
            # Если число Маха М > 1, то необходимо увеличить промежуточное давление
            elif M_i < 1:
                # если число Маха изменятся линейно
                if not(np.abs(previous_M_lower - M_i) < 0.1):
                    p_i = np.mean([p_i, self.__p_1])
                # если число Маха перестало линейно изменяться
                else:
                    p_i = 0.98 * p_i
                previous_M_lower = M_i
                
            # Если число Маха М < 1, то необходимо уменьшить промежуточное давление
            elif M_i > 1:                
                # если число Маха изменятся линейно
                if not(np.abs(previous_M_higher - M_i) < 0.1):
                    p_i = np.mean([self.__fp_0, p_i])
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
         
    def __get_normal_vane_area_and_alpha_1(self):
        """
            Функция предназначена для расчета площади сопловой решетки и угла выхода при числе Маха М < 1.
        """ 
        # Абсолютный угол выхода потока из сопловой решетки, в градусах
        alpha_1 = self.__alpha_1 = np.rad2deg(np.arcsin(self.__mu_1 * np.sin(np.deg2rad(self.__alpha_1eef))/self.__phi))
        
        # Выходная площадь сопловой решетки
        F_1 = self.__G_0 * self.__v_1t / self.__mu_1 / self.__c_1t
    
        return F_1, alpha_1
    
    def __get_critical_vane_area_and_alpha_1(self):
        """
            Функция предназначена для расчета площади сопловой решетки и угла выхода при числе Маха М > 1.
        """ 
        # Поиск критической точки в сопловой решетке
        Point_cr = self.__find_critical_point_2()

        # Критический удельный объем
        v_1cr = Point_cr.v(0)
                
        # Критическая абсолютная скорость в сопловой решетке
        c_1cr = np.sqrt(2000*(self.__fh_0 - Point_cr.h(0))) 
        
        # Приведенная относительная площадь
        q = self.__c_1t * v_1cr / c_1cr / self.__v_1t
                
        # Абсолютный угол выхода потока из сопловой решетки, в градусах
        alpha_1 = np.rad2deg((np.arcsin(self.__mu_1 * np.sin(np.deg2rad(self.__alpha_1eef)/self.__phi/q))))
        
        # Выходная площадь сопловой решетки
        F_1 = self.__G_0 * v_1cr / self.__mu_1 / c_1cr
        
        return F_1, alpha_1
        
    def __calculate_vane_length(self, phi, mu_1):
        '''
            Функция рассчитывает длину сопловых лопаток на основе коэффициентов скорости и расхода.
            
            Проблематика вычисления длины лопатки заключается в том, что необходимо определить длину
            без среднего диаметра. Для этого воспользуемся логичным выражением: 
            d = d_hub + Delra_pr + l_1. 
            
            Если подставить диаметр в выражение для определения высоты лопатки,
            то получится:
            l_1*e*(d_hub+Delta_pr+l_1) = F_1 / (pi * sin(alpha_1))
            
            Из чего следует, что необходимо решить квадратное уравнение. Для этого введем два коэффицинта
            b = e * d_hub + Delta_pr, c = F_1 / (pi * sin(alpha_1))
            
            Так уравнение примет привичный вид:
            l_1**2 + b*l_1 - c = 0
              
        '''
        
        self.__phi  = phi
        self.__mu_1 = mu_1
        
        # Коэффициент скорости для сопловой решетки
        # начальное приближение
            
        
        
        # Полный располагаемый теплоперепад ступени
        self.__fH_0 = self.__H_0 + np.power(self.__c_0, 2)/2000          
        
        # Полный располагаемый теплоперепад сопловой решетки
        self.__fH_0s = self.__fH_0*(1-self.__reaction)
        
        # Абсолютная теоретическая скорость на выходе из сопловой решетки
        self.__c_1t  = np.sqrt(2*self.__fH_0s*1000)
        
        # Термодинамическая точка f(p_0, t_0), статические параметры
        Point_0 = ThPoint(p=self.__p_0, t=self.__t_0)
        
        # Статическая величина энтальпии на входе в ступень
        h_0 = Point_0.h(0)
        
        # Статическая величина энтропии на входе в ступень
        s_0 = Point_0.s(0)
        
        
        ## Определение полных параметров на входе в сопловую решетку
        # Энтальпия полного торможения на входе в ступень
        self.__fh_0 = h_0 + np.power(self.__c_0, 2)/2000
        # Этропия полного торможения на входе в ступень
        self.__fs_0 = s_0
        # Термодинамическая точка полного торможения на входе в сопловую решетку
        fPoint_0 = ThPoint(h=self.__fh_0,s=self.__fs_0)
        # Давление полного торможения на входе в ступень
        self.__fp_0 = fPoint_0.p(0)
        
        # Определение тероеретических статических параметров на выходе из сопловой решетки
        # Теоретическая энтальпия на выходе из сопловой решетки
        h_1t = self.__fh_0 - self.__fH_0s
        # Термодинамическяа точка f=(h=h_1t, s=s_1t=s_0)
        Point_1t = ThPoint(h=h_1t, s=s_0)
        # Давление после сопловой решетки
        self.__p_1 = Point_1t.p(0)
        # Местная скорость звука на выходе из сопловой решетки 
        self.__a_1 = Point_1t.a(0)
        # Число маха для сопловой решетки
        self.__M_1t = self.__c_1t/self.__a_1
        # Удельный объем на выходе из сопловой решетки
        self.__v_1t = Point_1t.v(0)
        
        # Фиктивная скорость
        self.__c_f = np.sqrt(2000*self.__fH_0)
        print('c_f', self.__c_f)
            
        
        print('u', self.__u)
        
        # Выходная площадь сопловой решетки
        # Если число Маха M < 1
        if (self.__M_1t < 1):
            print("M < 1")
            # Выходная площадь сопловой решетки & Абсолютный угол выхода потока из сопловой решетки
            self.__F_1, self.__alpha_1 = self.__get_normal_vane_area_and_alpha_1()
                
        # Если число Маха М > 1
        else:
            print("M > 1")
            # Выходная площадь сопловой решетки & Абсолютный угол выхода потока из сопловой решетки
            self.__F_1, self.__alpha_1 = self.__get_critical_vane_area_and_alpha_1()
        
        # Если ступень рассчитывается на максимальный лопаточный КПД
        if (self.__is_stage_optimal):
            # Значение (u/c_ф)_опт
            self.__X = phi * np.cos(np.radians(self.__alpha_1eef))/(2*np.sqrt(1-self.__reaction))   
              
        # Если ступень рассчитывается НЕ на максимальный лопаточный КПД
        else:
            pass
            
        # Окружная скорость
        self.__u = self.__X * self.__c_f
        
        # Средний диаметр
        self.__d = self.__u / np.pi / self.__n
        print("d", self.__d)
            
        # Выходная высота сопловых лопаток 
        self.__l_1 = self.__F_1/(np.pi * self.__d * np.sin(np.deg2rad(self.__alpha_1eef)))
        print('l_1', self.__l_1)
        
        
        # Если лопатки получились меньше чем 0.013 мм
        if (self.__l_1 < 0.013):
            # Степень парциальности != 1
            self.__partiality = 6 * np.sqrt(self.__l_1)
            self.__l_1 = self.__l_1 / self.__partiality
        else:
            # Степень парциальности = 1
            self.__partiality = 1
            
    def __calculate_vane_speed_triangle(self):
        """
            Функция рассчитывает величины соплового треугольника скоростей.
            Список рассчитываемых величины:
            с_1         - Действительная абсолютная скорость на выходе из сопловой решетки
            beta_1      - Угол выхода из сопловой решетки в относительном движении
            w_1t,w_1    - Относительная теоретическая и действительная скорость на входе в рабочую решетку
        """
        # Действительная абсолютная скорость на выходе из сопловой решетки
        self.__c_1 = self.__c_1t * self.__phi
        
        # Действительная относительная скорость на входе в рабочую решетку
        self.__w_1 = np.sqrt(np.power(self.__c_1, 2) + np.power(self.__u, 2) - 2 * self.__u * self.__c_1 * np.cos(np.deg2rad(self.__alpha_1)))
        
        # Теоретическая относительная скорость на входе в рабочую решетку
        self.__w_1t = self.__w_1 / self.__psi
        
        # Угол выхода из сопловой решетки в относительном движении
        self.__beta_1 = self.__alpha_1 + np.rad2deg(np.arccos((np.power(self.__c_1, 2) + np.power(self.__w_1, 2) - np.power(self.__u, 2))/(2*self.__c_1*self.__w_1))) 
    
    def __calculate_vane_losses(self, phi):
        self.__Delta_Hs = np.power(self.__c_1t, 2) / 2000 * (1 - np.power(phi, 2))
    
    
    def __calculate_rotor_length(self):
        """
            Рассчет длины рабочей лопатки
        """ 
        self.__l_2 = self.__l_1 + self.__Delta_pr
    
    def __calculate_rotor_speed_triangle(self, psi, mu_2):
        self.__psi = psi
        self.__mu_2 = mu_2
        
        # Теоретиская энтальпия после сопловой решетки
        self.__h_1t = self.__fh_0 - self.__fH_0s
        
        # Статическая энтальпия на входе в рабочую решетку
        self.__h_1 = self.__h_1t + self.__Delta_Hs
        
        # Термодинамиская точка f(p=p_1, h_1)
        Point_1 = ThPoint(p=self.__p_1, h=self.__h_1)
        
        # Энтропия на входе в рабочую решетку
        s_1 = Point_1.s(0)
        fs_1 = s_1
        
        # Термодинамическая точка (точка соответствующая 2`t)
        Point_2tt = ThPoint(h=self.__fh_0 - self.__fH_0, s=self.__fs_0)
        
        # Статическое давление на выходе из ступени
        self.__p_2 = Point_2tt.p(0)
        
        # Термодинамическая точка на выходе из ступени (2t)
        Point_2t = ThPoint(p=self.__p_2, s=fs_1)
        
        # Теоретичская энтальпия на выходе из ступени
        self.__h_2t = Point_2t.h(0)
        
        # Располагаемый теплоперепад рабочей решетки
        H_0r = self.__h_1 - self.__h_2t
        
        # Полный располагаемый теплоперепад ступени
        self.__fH_0r = H_0r + np.power(self.__w_1t, 2)/2000
         
        # Теоретическая относительная скорость в рабочей решетке
        self.__w_2t = np.sqrt(2000*self.__fH_0r)
        
        # Действительная относительная скорость в рабочей решетке
        self.__w_2 = self.__w_2t * self.__psi
        
        
         
        
    def __calculate_rotor_losses(self, psi):
        self.__Delta_Hr = np.power(self.__w_2t, 2) / 2000 * (1 - np.power(psi, 2))
        
    
    def __get_empirical_phi_mu_1(b_1,l_1):
        """ Возвращает рассчитанный на основе эмпирических зависимостей 
            коэффициент скорости и расхода для сопловой решетки. """
        
        return (0.98 - 0.008*(b_1/l_1)),(0.982 - 0.005*(b_1/l_1))
    
    def __get_empirical_psi_mu_2(b_2,l_2):
        """ Возвращает рассчитанный на основе эмпирических зависимостей 
            коэффициент скорости и расхода для рабочей решетки. """
            
        return (0.96 - 0.014*(b_2/l_2)),(0.965-0.01*(b_2/l_2))
    
    
    
    # TODO
    def __check_bending_stresses(self):
        pass
    

    def __calculate_vane(self, phi, mu_1, b_1):
        '''
        Расчет сопловой решетки.
        
        Алгоритм состоит из следующих пунктов:
        1. Расчет длины сопловых лопаток
        2. Расчет 
        '''
        # Рассчитываем длину лопатки
        self.__calculate_vane_length(phi=phi, mu=mu_1)
        
        # Рассчитываем треугольник скоростей для сопловой решетки
        self.__calculate_vane_speed_triangle()
        
        # Рассчитываем потери в сопловой решетке
        self.__calculate_vane_losses(phi=phi)
    
        
        ## Определение давления после сопловой решетки
        # Степень реактивности на среднем диаметре
        #reaction = reaction_k + 1.8/
        

    def __calculate_rotor(self):
        '''
        Расчет рабочей решетки
        '''

    def __calculate_losses(self, c1, c2, alpha1, alpha2, l):
        '''
        Расчет потерь
        '''
        return 0.04
    
    def calculate(self):
        '''
        Проводим расчет ступени
        Примерный алгоритм расчета:
        1. Задаются хорды
        2. Задаются коэффициенты скорости и расхода для сопловой решетки
        3. Просчитывается сопловая решетка
        4. Просчитывается рабочая решетка
        5. Уточняются коэффициенты скорости и расхода, при не совпадении в три знака перед запятой
        алогритм начинается с пункта 2. но с уточненными коэффициентами скорости и расхода
        5. Просчитывается треугольник скоростей для сопловой решетки
        6. Просчитывается треугольник скоростей для рабочей решетки
        7. Подбирается профиль для рабочей решетки 
        8. Рабочая лопатка проверяется на изгиб, при невыполнении условий возвращается к пункту 1, хорда увеличивается (но как?)
        9. ...
        '''
        # Первое минимальное приближение хорд лопаток
        # 1.
        b_1 = 0.030
        b_2 = 0.020
        # 2.
        phi = 0.975
        mu_1 = 0.985
        
        while (True):
            # 3
            self.__calculate_vane(phi=phi, mu_1=mu_1, b_1=b_1)
            
            # 4
            self.__calculate_rotor(psi=psi, )
            
        
        
        

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
            'xi_s':     self.__xi_s, 
            'xi_r':     self.__xi_r, 
            'eff_ol':   self.__eff_ol,
            'E_0':      self.__E_0
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
        b:float,                              # [мм]      Хорда 
        B:float,                              # [мм]      Ширина
        min_sec_length:float,                 # [мм]      Горло
        rel_t:float,                          # []        Относительный шаг
        ins_angle:float,                      # [deg]     Установочный угол
        in_edge:float,                        # [мм]      Входная кромка
        out_edge:float,                       # [мм]      Выходная кромка
        #------------------------------------ # ---------------------------------------------------------------
        # Прочностные характеристики   
        # ----------------------------------- # ---------------------------------------------------------------
        x_weight_center:float,                # [мм]      x-координата центра тяжести
        y_weight_center:float,                # [мм]      y-координата центра тяжестт
        profile_area:float,                   # [см^2]    Площадь профиля
        inertia_moment_xx:float,              # [см^4]    Момент инерции относительно ХХ
        inertia_moment_yy:float,              # [см^4]    Момент инерции относительно YY
        resistance_moment_xx_back:float,      # [см^3]    Момент сопротивления относительно ХХ, по спинке
        resistance_moment_xx_edge:float,      # [см^3]    Момент сопротивления относительно ХХ, по кромкам
        resistance_moment_yy_in_edge:float,   # [см^3]    Момент сопротивления относительно YY, по вх. кромке
        resistance_moment_yy_out_edge:float,  # [см^3]    Момент сопротивления относительно YY, по вых. кромке
        
        #ins_angle_exp,
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
    

class GridProfile:
    """
        Класс, храняющий основные параметры турбинных сопловых и рабочих решеток МЭИ.
    """
    __data = pd.read_excel('API/grids.xlsx')
    
    def __init__(self, name:str):
        try:
            # Копируем необходимую строку параметров для выбранной решетки
            profile_parameters = self.__data[self.__data['Name'] == name]
        
            # Если строка найдена не одна или не найдена
            if (len(profile_parameters) != 1):
                raise Exception("Invalid number of founded grids, number of greeds: {}".format(len(profile_parameters)))

            self.__name = name
            self.__existing_grid = ExistingGrid(
                b                               = float(profile_parameters['b']),
                B                               = float(profile_parameters['B']),
                min_sec_length                  = float(profile_parameters['a']),
                rel_t                           = float(profile_parameters['rel_t']),
                ins_angle                       = float(profile_parameters['ins_angle']),
                in_edge                         = float(profile_parameters['Delta_in']),
                out_edge                        = float(profile_parameters['Delta_out']),
                x_weight_center                 = float(profile_parameters['center_X']),
                y_weight_center                 = float(profile_parameters['center_Y']),
                profile_area                    = float(profile_parameters['F']),
                inertia_moment_xx               = float(profile_parameters['Ixx']),
                inertia_moment_yy               = float(profile_parameters['Iyy']),
                resistance_moment_xx_back       = float(profile_parameters['Wxx_back']),
                resistance_moment_xx_edge       = float(profile_parameters['Wxx_edge']),
                resistance_moment_yy_in_edge    = float(profile_parameters['Wyy_in_edge']),
                resistance_moment_yy_out_edge   = float(profile_parameters['Wyy_out_edge'])
            )
            
        except Exception as ex:
            print(*ex.args)
            
        # Getters      
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
        
        
        
        
    
    
