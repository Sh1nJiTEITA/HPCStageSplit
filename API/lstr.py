# Далее идет класс для вывода в LaTeX'е
import pint
import numpy as np

class lstr(str):
    
    # классическая версия str.format с изменением замены {} на [] и без поддержки [1] [2]
    def format(self, *args, **kwargs):
        
        for i in args:
            self = self.replace("[]", lstr(i), 1)
        
        for i in kwargs.items():
            self = self.replace(lstr('[' + i[0] + ']'), lstr(i[1]))
        
        return lstr(self)
    
    def format_comment(self, *args, **kwargs):
        
        for i in args:
            self = self.replace("%[]", lstr(i), 1)
        
        for i in kwargs.items():
            self = self.replace(lstr('%[' + i[0] + ']'), lstr(i[1]))
        
        return lstr(self)
    
    
    def __correct_round(self, value, ALIGN_MODE = 0):
        '''
            Правильное округление величин, достаточное для инженерных расчетов:
            
            # if ==0,1 : 0.0044, .0044 -> 0.004
            # if ==2: 12.3313123 -> 12.33
            # if >=3: 343.131233 -> 343.1, 34141414.54 ->34141414.5
            # if 0,1 : 0.00, .0
            
        '''    
        
        dot_pos = lstr(value).find(".")
        
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
            for i in lstr(value)[dot_pos + 1:]:
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
            Версия lstr.format(...) с поддержкой единиц измерения с использоваением округления
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

        
    def dformat(self, inlist=[],*args, **kwargs):
        '''
            Версия lstr.format(...) с поддержкой единиц измерения с использоваением округления
        '''
        
        _list = list()
        for i in inlist:
            if (type(i) == np.ndarray):
                raise Exception("this:{} is np.ndarray")
            if (isinstance(i[1], pint.Quantity)):
                if (type(i[1].m) == np.ndarray):
                    continue
                _list.append((i[0], self.__correct_round(i[1].m, ALIGN_MODE=0)))
            else:
                _list.append((i[0], i[1]))  
        
        _args = list()
        for i in args:
            if (type(i) == np.ndarray):
                raise Exception("this:{} is np.ndarray")
            if (isinstance(i, pint.Quantity)):
                if (type(i.m) == np.ndarray):
                    continue
                _args.append(self.__correct_round(i.m, ALIGN_MODE=0))
            else:
                _args.append(i)
                
        _kwargs = list()
        for i in kwargs.items():
            # if (type(i[1]) == np.ndarray):
            #     raise Exception("this:{} is np.ndarray")
            if (isinstance(i[1], pint.Quantity)):
                if (type(i[1].m) == np.ndarray):
                    continue
                _kwargs.append((i[0], self.__correct_round(i[1].m, ALIGN_MODE=0)))
            else:
                _kwargs.append((i[0], i[1]))  
                
        return self.format(*tuple(_args), **dict(_kwargs), **dict(_list))
    
    def get_large(self):
        local_lstr = self
        local_lstr = lstr(r"\Large{") + local_lstr + lstr(r"} \\ \ \\")
        return local_lstr

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

a = lstr("[][][b][c]")

#print(a.format(1,2,b=" 13",c=" vd"))

#print("{}{}{b}".format("AAA", "b", b="3"))
#swap_br_in_math_string("3213123", a=10, b=20)
#Math(lstr("\Delta [H_0] = [H_0v]\ кДж/кг").format(H_0="H_0",H_0v=33))

#print(type(lstr(r"3.2 \ \theta = []").dformat(theta)))