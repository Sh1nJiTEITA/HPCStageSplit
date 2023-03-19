# Используемые библиотеки

# Общие
import numpy as np

# Собственные
from . import ThPoint

# Константы
def DRY_STEAM_CRITICAL_PRESSER_RATIO(): return 0.5774
def OVERHEATED_STEAM_CRITICAL_PRESSER_RATIO(): return 0.5447



class Stage:
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
        
        self.__is_stage_optimal:bool = isStageOptimal
    
        
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
        self.__X        = None          # []         Значение (u/c_f)_отп
        self.__Delta_Hs = None          # [кДж/кг]   Потери в сопловой решетке
        self.__Delta_Hr = None          # [кДж/кг]   Потери в рабочей решетке
        self.__xi_s     = None          # []         Относительные потери в сопловой решетке
        self.__xi_r     = None          # []         Относительные потери в рабочей решетке
        self.__eff_ol   = None          # []         Относительный лопаточный КПД ступени
        self.__E_0      = None          # [кДж/кг]   Располагаемая энергия ступени
        
        # Промежуточные величины
        #
        
        self.__fH_0     = None          # [кДж/кг]   Полный располагаемый теплоперепад ступени
        self.__fH_0s    = None          # [кДж/кг]   Полный располагаемый теплоперепад сопловой решетки
        self.__c_1t     = None          # [м/с]      Абсолютная теоретическая скорость на выходе из сопловой решетки
        self.__fH_0r    = None          # [кДж/кг]   Полный располагаемый теплоперепад Рабочей решетки
        self.__p_1      = None          # [МПа]      Давление на выходе из сопловой решетки
        self.__M_1t     = None          # []         Число Маха на выходе из сопловой решетке (в горле)
        self.__mu_1     = None          # []         Коэффициент расхода для сопловой решетки
        self.__fp_0     = None          # [МПа]      Давление полного торможения на входе в ступень
        self.__fh_0     = None          # [кДж/кг]   Энтальпия полного торможения на входе в ступень
        self.__fs_0     = None          # [кДж/кг/К] Энтропия полного торможения на входе в ступень
        self.__alpha_1  = None          # [deg]      Абсолютный угол выхода потока из сопловой решетки
        self.__partiality = None        # []         Степень парциальности
        self.M_i       = None
         
        self.__define_speed_koef()
    
    
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
         
        
        
    
    def __define_speed_koef(self):
        '''
            Функция задает и уточняет коэффициенты скорости.
            Расчет отличается в зависимости от типа ступени: рассчитывается ли ступень на максимальный КПД, или нет
        '''
        
        # Если ступень рассчитывается на максимальный лопаточный КПД
        if (self.__is_stage_optimal):
            
            # Коэффициент скорости для сопловой решетки
            # начальное приближение
            phi = 0.96
            
            # Необходимо определить вид течение
            
            
            # Значение (u/c_ф)_опт
            X = phi * np.cos(np.radians(self.__alpha_1eef))/(2*np.sqrt(1-self.__reaction))
            
            # Полный располагаемый теплоперепад ступени
            self.__fH_0 = self.__H_0 + np.power(self.__c_0, 2)/2000
            
            # !Средний диаметр ступени
            
            # Полный располагаемый теплоперепад сопловой решетки
            self.__fH_0s = self.__fH_0*(1-self.__reaction)
            
            # !Полный располагаемый теплоперепад рабочей решетки
            self.__fH_0r = self.__fH_0 - self.__fH_0s
            
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
            v_1t = Point_1t.v(0)
            
            # !!! Задаемся коэффициентом расхода в сопловой решетке
            self.__mu_1 = 0.99
            
            
            
            # Выходная площадь сопловой решетки
            # Если число Маха M < 1
            if (self.__M_1t < 1):
                print("M < 1")
                # Выходная площадь сопловой решетки
                F_1 = self.__G_0 * v_1t / self.__mu_1 / self.__c_1t
                
                 # Абсолютный угол выхода потока из сопловой решетки
                self.__alpha_1 = np.arcsin(self.__mu_1 * np.sin(np.deg2rad(self.__alpha_1eef))/phi)
                
            # Если число Маха М > 1
            else:
                print("M > 1")
                # Критическая термодинамическая точка
                Point_cr = self.__find_critical_point_2()
                
                # ----------------------------------------------
                print(Point_cr)
                print("eps_cr_div", Point_cr.p(0)/self.__fp_0)
                k = Point_cr.k(0)
                eps_kr = np.power(2/(k+1),k/(k-1))
                print("eps_cr",eps_kr)
                # ----------------------------------------------
                
                # Критический удельный объем
                v_1cr = Point_cr.v(0)
                
                # Критическая абсолютная скорость в сопловой решетке
                c_1cr = np.sqrt(2000*(self.__fh_0 - Point_cr.h(0))) 
                
                # Выходная площадь сопловой решетки
                F_1 = self.__G_0 * c_1cr / self.__mu_1 / c_1cr
                
                # Приведенный расход в критическом сечении
                q = self.__c_1t * v_1cr / c_1cr / v_1t
                print("q", q)
                
                # Абсолютный угол выхода потока из сопловой решетки, в градусах
                self.__alpha_1 = np.rad2deg((np.arcsin(self.__mu_1 * np.sin(np.deg2rad(self.__alpha_1eef)/phi/q))))
                print('alpha_1', self.__alpha_1)
                
                
            # Выходная высота сопловых лопаток 
            self.__l_1 = F_1/(np.pi * (self.__d_hub + self.__Delta_pr + self.__l_1) * np.sin(np.deg2rad(self.__alpha_1)))
            
            # Если лопатки получились меньше чем 0.013 мм
            if (self.__l_1 < 0.013):
                # Степень парциальности != 1
                self.__partiality = 6 * np.sqrt(self.__l_1)
                self.__l_1 = self.__l_1 / self.__partiality
            else:
                # Степень парциальности = 1
                self.__partiality = 1    
              
            
            
        # Если ступень рассчитывается НЕ на максимальный лопаточный КПД
        else:
            pass
            
        
    

    def __calculate_vane(self, phi):
        '''
        Расчет сопловой решетки
        '''
        
        ## Статические параметры на входе в ступень
        # Термодинамическая точка f(p_0, t_0), статические параметры
        Point_0 = ThPoint(p=self.__p_0, t=self.__t_0)
        
        # Статическая величина энтальпии на входе в ступень
        h_0 = Point_0.h(0)
        
        # Статическая величина энтропии на входе в ступень
        s_0 = Point_0.h(0)
        

        ## Полные параметры на входе в ступень
        # Полная величина энтропии на входе в ступень
        fs_0 = s_0
        
        # Полная величина энтальпи на входе в ступень
        fh_0 = h_0 + np.power(self.__c_0, 2)/2000
        
        # !Термодинамическая точка f(fp_0, fs_0), полные параметры
        fPoint_0 = ThPoint(h=fh_0,s=fs_0)
        
        
        ## Определени давления после рабочей решетки
        # Статическая теоретическая величина энтальпии после рабочей решетки
        h_2t = h_0 + self.__H_0
        
        # Термодинамическая точка f(h_2t, s_0), статические параметры
        Point_2 = ThPoint(h=h_2t, s=s_0)
        
        # !Давление на выходе из рабочей решетки
        p_2 = Point_2.p(0)
        
        ## Определение 
        
        
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
        '''

    def get_results(self):
        return {
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