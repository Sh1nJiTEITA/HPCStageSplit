import pint
import numpy as np

from . import ThPoint
# from . import Stage

class tlstr(str):
    
    # классическая версия str.format с изменением замены {} на [] и без поддержки [1] [2]
    def format(self, *args, **kwargs):
        
        for i in args:
            self = self.replace("<<>>", tlstr(i), 1)
        
        for i in kwargs.items():
            self = self.replace(tlstr('<<' + i[0] + '>>'), tlstr(i[1]))
        
        return tlstr(self)
    
    def format_comment(self, *args, **kwargs):
        
        for i in args:
            self = self.replace("%<<>>", tlstr(i), 1)
        
        for i in kwargs.items():
            self = self.replace(tlstr('%<<' + i[0] + '>>'), tlstr(i[1]))
        
        return tlstr(self)
    
    
    def __correct_round(self, value, ALIGN_MODE = 0):
        '''
            Правильное округление величин, достаточное для инженерных расчетов:
            
            # if ==0,1 : 0.0044, .0044 -> 0.004
            # if ==2: 12.3313123 -> 12.33
            # if >=3: 343.131233 -> 343.1, 34141414.54 ->34141414.5
            # if 0,1 : 0.00, .0
            
        '''    
        if not(isinstance(value, (float, int))):
            return value
        
        elif isinstance(value, ThPoint):
            return value
        elif value is None:
            return value
        elif isinstance(value, str):
            return value
        
        
        dot_pos = tlstr(value).find(".")
        
        if 'e' in str(value):
            exp_mode = True
            exp_degree = str(value)[str(value).find('e') + 1:]
            value = float(str(value)[:str(value).find('e')])
            
            for pos,letter in enumerate(exp_degree):
                if letter == '-':
                    continue
                elif letter == '+':
                    continue
                elif letter == '0':
                    continue
                else:
                    main_number_pos = pos
                    break
            exp_degree = '\cdot 10^{'+ exp_degree[0] + exp_degree[pos:] + '}'
        else:
            exp_mode = False
            exp_degree = ''        
        
        
        if ((dot_pos == 0) or (dot_pos == 1)):
            # Если включен режим округления для ровного вырванивния
            # Так как число может начинаться с любого количества нулей, а нам нужно только последние три цифры после, то
                # необходимо найти координаты начала этих значимых цифр
            last_null_pos = 0
            
                # ищем
            for i in tlstr(value)[dot_pos + 1:]:
                if (i != '0'):
                    break
                else:
                    last_null_pos += 1
            
            if (ALIGN_MODE):
                if (len(str(value)) < 2 + last_null_pos+3):
                    return str(value) + exp_degree
            
                # округляем до позиции первой цифры и еще две цифры сверху (в сумме 3)
                out_value = round(value, last_null_pos+3)
            
                # Если число округляется так, что значащих цифр станвится 2 (и еще ноль, который не отрисовывается), то 
                # добавляем его для ровной отрисовки далее
                # 0.000010999 -> 0.0000110 а не 0.000011
            
                # литералы ниже: 2-это "0." в начале числа, 3-это последние значащие цифры
                if len(str(out_value)) != 2 + last_null_pos+3:
                    return str(out_value) + '0'+ exp_degree
                else:
                    return str(out_value) + exp_degree
            # режим выравнивания выключен
            else:
                return str(round(value, last_null_pos+3)) + exp_degree
            
        elif (dot_pos == 2):
            if (ALIGN_MODE):
                out_value = round(value, 2)
            
                # Если число не состоит из двух первых знаков, запятой и двух последующих, то добавляем ноль для ровной отрисовки
                # 14.0999 -> 14.10 а не 14.1
            
                # литералы ниже: 2 - первые два знака до запятой, 1 - точка между целыми и дробными, 2 - два знака после запятой
                if len(str(out_value)) != 2 + 1 + 2:
                    return str(out_value) + '0' + exp_degree
                else:
                    return str(out_value) + exp_degree
            else:
                return str(round(value,2)) + exp_degree
            
        elif (dot_pos >= 3):
            if (ALIGN_MODE):
                out_value = round(value, 1)
            
                # Если большое число округлилось в большую сторону, то добавляем ноль для правильной отрисовки далее
                # 19999.9999 -> 20000.0 а не 20000
            
                # если длина строки меньше, чем была бы с знаком после точки
                if len(str(out_value)) <= dot_pos+1:
                    return str(out_value) + '0' + exp_degree
                else:
                    return str(out_value) + exp_degree
            else:
                return str(round(value,1)) + exp_degree
        else:
            return str(value) + exp_degree
    
    def dlformat(self, *args, **kwargs):
        '''
            Версия tlstr.format(...) с поддержкой единиц измерения с использоваением округления
        '''
        
        _args = list()
        for i in args: 
            if (isinstance(i, pint.Quantity)):
                _args.append(self.__correct_round(i.m, ALIGN_MODE=1))
            else:
                _args.append(i)
                
        _kwargs = list()
        for i in kwargs.items():
            if (isinstance(i[1], pint.Quantity)):
                _kwargs.append((i[0], self.__correct_round(i[1].m, ALIGN_MODE=1)))
            else:
                _kwargs.append((i[0], i[1]))  
                
        return self.format(*tuple(_args), **dict(_kwargs))
    
    def performat(self, *args, **kwargs):
        _args = list()
        for i in args: 
            if (isinstance(i, pint.Quantity)):
                _args.append(self.__correct_round(i.m, ALIGN_MODE=0))
            else:
                _args.append(i)
                
        _kwargs = list()
        for i in kwargs.items():
            if (isinstance(i[1], pint.Quantity)):
                _kwargs.append((i[0], self.__correct_round(i[1].m, ALIGN_MODE=0)))
            else:
                _kwargs.append((i[0], i[1]))  
                
        return self.format_comment(*tuple(_args), **dict(_kwargs))

        
    def dformat(self, inlist=[],mode=0,*args, **kwargs):
        '''
            Версия tlstr.format(...) с поддержкой единиц измерения с использоваением округления
        '''
        
        _list = list()
        for i in inlist:
            if (type(i[1]) == np.ndarray):
                continue
            if (type(i[1]) == ThPoint):
                continue
            if (type(i[1]) == bool):
                continue
            # if (type(i[1]) == None):
            #     continue
            # if (type(i[1]) == str):
            #     _list.append((i[0], i[1]))
            #     continue
            if (isinstance(i[1], pint.Quantity)):
                if (type(i[1].m) == np.ndarray):
                    continue
                _list.append((i[0], self.__correct_round(i[1].m, ALIGN_MODE=0)))
            else:
                if not mode:
                    _list.append((i[0], i[1]))
                else:
                    _list.append((i[0], self.__correct_round(i[1], ALIGN_MODE=0)))
        
        _args = list()
        for i in args:
            if (type(i) == np.ndarray):
                continue
            if (isinstance(i, pint.Quantity)):
                _args.append(self.__correct_round(i.m, ALIGN_MODE=0))
            else:
                _args.append(i)
                
        _kwargs = list()
        for i in kwargs.items():
            
            if (isinstance(i[1], pint.Quantity)):
                if (type(i[1].m) == np.ndarray):
                    continue
                _kwargs.append((i[0], self.__correct_round(i[1].m, ALIGN_MODE=0)))
            else:
                if not mode:
                    _list.append((i[0], i[1]))
                else:
                    _list.append((i[0], self.__correct_round(i[1], ALIGN_MODE=0)))  
                
        return self.format(*tuple(_args), **dict(_kwargs), **dict(_list))
    
    def get_large(self):
        local_tlstr = self
        local_tlstr = tlstr(r"\Large{") + local_tlstr + tlstr(r"} \\ \ \\")
        return local_tlstr

def get_dict_of_class(input_class):
    vec = {}
    for item in input_class.__dict__:
        key:str = str(item)
        #print(input_class.__class__.__name__)
        key = key.replace(input_class.__class__.__name__, '')
        key = key.replace('_', '')
        #print(key)
        
        vec[key] = input_class.__dict__[item]
        
    return vec

a = tlstr("[][][b][c]")