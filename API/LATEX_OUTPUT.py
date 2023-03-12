from IPython.display import Latex, display
import numpy as np
from . import TCPv2
from pint import get_application_registry, Quantity
u = get_application_registry()
Q = Quantity


def make_math_string(input_str:str, **args):
    for i in args:
        print(i)




def output_value_var_str(var, value):
    if not((isinstance(value, int))or
           (isinstance(value, float))or
           (isinstance(value,np.float64))or
           (isinstance(value,np.int32))):
        value = "{}\ {}".format(value.magnitude, value.units)
    #return (f"""\\begin{{align}}{var} = {value}\\end{{align}}""")
    return (f"""{var} = {value}""")
    #  return Latex(f"""\\begin{{align}}
                 #\\{var} = {value}\\end{{align}}""")
def op_latex_output(data):
    begin_str = f'\\begin{{dcases}}'
    end_str = f'\\end{{dcases}}'
    output_str = ""
    
    if not(isinstance(data, TCPv2.ThPoint)):
        for i in data:
            output_str+=(output_value_var_str(i[0], i[1]))
            output_str+=r" \\ "
    else:
        #print("y")
        param_container = data.dl()
        for i in param_container:
            output_str+=(output_value_var_str(i[0], i[1]))
            output_str+=r" \\ "
        
    
    re_str = begin_str + output_str
    re_str += end_str
    
    return re_str


# TODO ПЕРЕПИСАТЬ !!!!!!!!!!!! (повтор кода, сделать классом)
def op_latex_output_2(data):
    '''
        Улучшенная версия (вторая) генерации вывода параметров в LaTeX
        
        На вход принимает массив и только массив. 
        1. Если значение одно - необходимо обложить его
        следующим образом: 
        
        op_latex_output_2([
            ("VARIABLE_NAME", value),
            ...
        ])
        
        2. Если входное значение является термодинамической точкой (TCPv2.TCPv2.ThPoint или его
        последующие версии), то в качестве VARIABLE_NAME используется индекс переменной, связанной
        с термодинамической точкой. В качестве значения на вход подается сама переменная с 
        разрешением TCPv2.TCPv2.ThPoint, то есть:
        
        Point_1t = TCPv2.TCPv2.ThPoint(p=10, t=550)
        op_latex_output_2([
            ("1t", Point_1t),
            ...
        ])
        
        3. Имеется возможность взаимной подачи переменных разных типов (первый случай) вместе с
        переменными типа TCPv2.TCPv2.ThPoint в одном листе. В этом случае генерация текста будет происходить
        в порядке перечисления в листе. Переменные термодинамической точки всегда генерируются в 
        одном и том же порядке подряд.
        То есть:
        
        Point_1t = TCPv2.TCPv2.ThPoint(p=10, t=550)
        op_latex_output_2([
            ("1t", Point_1t),
            ("VARIABLE_NAME", value),
            ...
        ])
        
        TODO: доделать более точное объявление.
        
        4. Имеется возможность вывода массивов 
        
        ADDITIONAL. Так же имеет возможно обозначения полных параметров (полного торможения) для 
        термодинамических точек (типа TCPv2.TCPv2.ThPoint), для этого необходимо добавить ключевое слово 'FULL'
        к индексу, в любое место. Очевидно, что невозможно использовать индекс с значением 'FULL'. 
        Обозначение полных параметров - черта над переменой (индекс не надчеркивается).
        То есть:
        
        Point_1t = TCPv2.TCPv2.ThPoint(p=10, t=550)
        op_latex_output_2([
            ("FULL1t", Point_1t),
            ("FULL 1t", Point_1t),
            ("1tFULL", Point_1t),
            ...
            
        Для полного погружения необходимо посмотреть функцию TCPv2.TCPv2.ThPoint.dl() и name_by_index_tex: dict.
        ])
    '''

    # Начало и конец строки в LaTeX, используется блок типа dcases (она же фигурная скобка)
    # такая же как, что и в системах
    begin_str = f'\\begin{{dcases}}'
    end_str = f'\\end{{dcases}}'
    
    # Строка, которая после отправится на выход из функции
    output_str = ""
    
    # Начинаем перебирать каждый элемент поданного массива
    # i = data[it] = ("VARIABLE_NAME", value)
    
    def find_br_pos(input_str):
        left_br_pos = input_str.rfind('_{') + 1
        right_br_pos = None
        br_counter = 1
        
        for it in range(left_br_pos, len(input_str)):
            if (input_str[it] == "}"):
                br_counter += 1
            else:
                br_counter -= 1
            if (br_counter == 0):
                right_br_pos = it
                break
        return (left_br_pos, right_br_pos)
    
    # def is_input_name_contains_origin(input_str):
    #     if (('ORIGIN'in input_str)):
    #         return True
    #     else:
    #         return False
    
    
        
        
    
    for i in data:
        # Далее алгоритм обрабывает обычные значения и термодинамические точки по-разному
        
        
        # Если итератор указывает на обычное значение (не TCPv2.ThPoint)
        if (not(isinstance(i[1], TCPv2.ThPoint))):
            ORIGIN = 0
            if ((isinstance(i[1], np.ndarray))):
                if ('ORIGIN' in i[0]):
                    ORIGIN = 1
               
                    
                if ('_{' in i[0]):
                    left_br_pos, right_br_pos = find_br_pos(i[0])
                    for j in range(0, len(i[1])):
                        local_name = i[0][:right_br_pos] + str(j + ORIGIN) + i[0][right_br_pos:]
                        output_str+=(output_value_var_str("{" + local_name + "}", i[1][j]) + r' \\ ')
                else:
                    for j in range(0, len(i[1])):
                        output_str+=(output_value_var_str("{" + i[0] + "}_{" + str(j + ORIGIN) + "}", i[1][j]) + r' \\ ')
                continue
                
            else:
                if (isinstance(i[1].m, np.ndarray)):
                    if ('ORIGIN' in i[0]):
                        ORIGIN = 1
                      
                    
                    if ('_{' in i[0]):
                        left_br_pos, right_br_pos = find_br_pos(i[0])
                        for j in range(0, len(i[1])):
                            local_name = i[0][:right_br_pos] + str(j) + i[0][right_br_pos:]
                            output_str+=(output_value_var_str("{" + local_name + "}", i[1][j]) + r' \\ ')
                    else:
                        for j in range(0, len(i[1])):
                            output_str+=(output_value_var_str("{" + i[0] + "}_{" + str(j + ORIGIN) + "}", i[1][j]) + r' \\ ')
                    
                    
            
                else:
                # Берется значение переменной VARIABLE_NAME, оборачивается в скобки
                # и отправляется в сборочную функцию (объявлена выше)
                    output_str+=(output_value_var_str("{" + i[0] + "}", i[1]))
                # Добавляется перенос строки, аналог '\n' для LaTeX
                    output_str+=r" \\ "
        
        # Если итератор указывает на значения типа TCPv2.ThPoint
        else:
            # Переменная, которая будет содержать текст термодинамической точки
            thpoint_str = ""
            
            # Переменная отвечающая за режим обозначения полных параметров
            # Если значение переменной False - полные параметры отключены
            # Если значение переменной True - полные параметры включены
            full_parameter_mode = False
            
            # Обновление переменной для обозначения полных параметров
            # Если есть ключевое слово, то переменная True, копируем значение индекса без ключевого слова
            if ('FULL' in i[0]):
                thpoint_index = i[0].replace('FULL','')
                full_parameter_mode = True
            # Если нет ключевого слова, то переменная False, копируем значение индекса без ключевого слова
            else:
                thpoint_index = i[0]
                full_parameter_mode = False
            
            # С помощью специальной функции из TCPv2.TCPv2.ThPoint получаем лист, состоящий из 
            # туплов (совокупность двух переменных), аналогичных ("VARIABLE_NAME", value),
            # к примеру:
            # param_container = i[1].dl() = [
            #     ("{p}",pVALUE),
            #     ("{h}",hVALUE),   
            # ]
            param_container = i[1].dl()
            
            # Далее начинаем перебирать данных контейнер для формирования строки параметров
            # термодинамических параметров
            # Не стоит путать этот перебор массива с перебором выше, который перебирал входные данные
            # Напоминаю, что входные данные могут состоять из нескольких обычных значений и точек типа TCPv2.ThPoint
            for j in param_container:

                # Дальнешая работа с переменными делится на два типа, которые в свою очередь еще на два. 
                
                # Первое разделение: Имеется ли уже индекс у параметра
                # Так как некоторые термодинамические параметры состоя не только из одной буквы,
                # например, C_p, C_v, k_v, для них необходимо "добавить" индекс, а не создать его, как для остальных
                
                # Второе разделение: На обозначение полных и неполных параметров
                # Для первых необходимо добавить кодовые значения в стиле LaTeX, чтобы обозначит черту сверху
                
                # Если переменная содержит '_' значит она нуждается в добавлении индекса, а не создании
                if ('_' in j[0]):
                    # Является ли входная точка полной (проверяется ранее)?
                    # Не является
                    if not(full_parameter_mode):
                        thpoint_str+=(output_value_var_str(j[0] + thpoint_index + '}', j[1]))
                    # Является
                    else:
                        thpoint_str+=(output_value_var_str("\=" + j[0] + thpoint_index + '}', j[1]))
                
                # Так как переменная не содержит '_', то для нее необходимо добавить индекс с нуля
                else:
                    # Является ли входная точка полной (проверяется ранее)?
                    # Не является
                    if not(full_parameter_mode):
                        thpoint_str+=(output_value_var_str(j[0] +"_{" + thpoint_index + "}", j[1]))
                    # Является
                    else:
                        thpoint_str+=(output_value_var_str("\=" + j[0] +"_{"+ thpoint_index + "}", j[1]))
                # После каждого параметр добавляем перенос строки
                thpoint_str+=r" \\ "
            
            # Добавляем все параметры к общей строке
            output_str += thpoint_str
            
            # Так как параметров много и надо их отделить, то добавляем перенос строки
            output_str += " \\ "
    
    # Форируем выходную строку из начальных и конечных строк и содержимого, сгенерированного в циклах
    # то есть совокупность термодинамических параметров точек и простых значений
    re_str = begin_str + output_str
    re_str += end_str
    
    re_str = re_str.replace('ORIGIN', '')
    
    return re_str


# функции для быстрого вывода в LaTeX для ipynb
def display_latex_ex_(data):
    display(Latex(op_latex_output(data)))

def display_latex_ex_2(data):
    display(Latex(op_latex_output_2(data)))
# A = TCPv2.ThPoint(p=10, t=550)
# print(op_latex_output_2(
#     [
#         ("FULL1f", A),
#         ("D", 550)
#     ]
# ))
# '''

# TEXT HEAR

# '''
# abc = Q(np.array([10, 30, 20, 40]), 'meter')
# print(abc)

# print(op_latex_output_2([
#     ('F', np.array([10,20,30])),
#     ('D', abc),
#     ("ddd", Q(10, "Hz")),
#     ("d1", TCPv2.ThPoint(p=4, t=550))
# ]))

# print(op_latex_output_2([
#     ('D ORIGIN', abc),
# ]))