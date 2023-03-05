import pint
import sympy as sm

import sympy.physics.units as un
import sympy.parsing.sympy_parser as sp


u = pint.get_application_registry()

# d_k = 10    * u("meter")
# l_21 = 20   * u("meter")
# d_21 = 0.9  * u("meter")
# v_2z = 33   * u('meter**3/kg')
# v_2t = 28   * u('meter**3/kg')

# s1 = pint.Quantity(33, "kJ/kg/K")
# s2 = pint.Quantity(1033, "kJ/kg/K")
def make_correct_str(data):
    return str(data.magnitude) + " * " + str(data.units)



def make_correct_str_base(data):
    '''
        Функция, получающая на вход объект (pint.QUantity, int, float) и преобразующее в необходимый
        для вычислений с учетом единиц измерения вид для использования совместно с Sympy.
        Для безопасности вычисления еденицы измерения приводятся к базовым.
        
        Пример: 
        
        1. a = pint.Quantity(99, 'j/kg')
        2. print(make_correct_str_base(a)) # prints: 99 * m**2/(c**2)/kg
        
    '''
    if not(isinstance(data,pint.Quantity)):
       return str(data)
    
    data = data.to_base_units()
    return str(data.magnitude) + " * (" + str(data.units) + ")"

def make_correct_str_base_(data):
    '''
        Функция, получающая на вход объект (pint.QUantity, int, float) и преобразующее в необходимый
        для вычислений с учетом единиц измерения вид для использования совместно с Sympy.
        Для безопасности вычисления еденицы измерения приводятся к базовым.
        
        Пример: 
        
        1. a = pint.Quantity(99, 'j/kg')
        2. print(make_correct_str_base(a)) # prints: 99 * m**2/(c**2)/kg
        
    '''
    if not(isinstance(data,pint.Quantity)):
       return str(data)
    
    data = data.to_base_units()
    return "(" + str(data.magnitude) + " * (" + str(data.units) + "))"

def swap_sqrt(data:str):
    '''
        Функция, получающая на вход строку (любое математическое выражение) преобразующее 
        все корни (sqrt(...)) в степенной вид (**0.5). 
        
        Пример. 
        Вызов функции с следующим аргументом: swap_sqrt('sqrt(3) / sqrt(4)').
        Вернет: '(3)**0.5 / (4)**0.5'
        
        Сложность данной задачи заключается в корректном трекинге скобок. Так как существуют
        усложненные выражения типа: 'sqrt(sqrt(sqrt())) - sqrt()' и похожее. Так как скробок 
        может быть очень много, необходимо неати именно ту, что отвечает за корень и добавить
        после нее степень.
        
        Найденное решение не является эффективным.
    '''
    
    # цикл while необходим для преобразования выражений, содежащих несколько корней
    # для каждого 'sqrt' будет проходит одна итерация цикла while
    while (True):
        # проверяем, содержится ли в входном выражении корни
        if ('sqrt' in data):
            
            # находим первый корень (поиск ведется до первого встреченного значения слева строки)
            pos = data.find('sqrt')
            
            # координата левой скобки, связанной с корнем
            # литера 4 - есть длина слово sqrt
            begin_br_pos = pos + 4

            # объявляем два вектора, которые будут содержать координаты (позиции)
            # левых и правых скобок в входной строке, что важно - до изменения
            left_br_poses = []
            right_br_poses = []
    
            # просматриваем каждый символ строки и записываем в вектора координаты скобок (всех)
            for i in range(0, len(data)):
                if (data[i] == '('):
                    left_br_poses.append(i)
        
                if (data[i] == ')'):
                    right_br_poses.append(i)

            # далее необходимо определить правую скобку, связанную с рассматриваемым
            # в итерации корнем
            # для этого создается вектор, содержащий позиции всех скобок, найденных ранее
            br_vec = []
            for i in range(0, len(left_br_poses)):
                br_vec.append(left_br_poses[i])
                br_vec.append(right_br_poses[i])
            
            # вектор сортируется так, чтобы позиции выставились по возрастанию
            br_vec.sort()
            
            # далее образуется еще один вектор, подражающий прошлому, но в каждой ячейке
            # хранится не только позиция, но и скобка (правая и левая)
            br_vec_ = []
            for i in range(0, len(br_vec)):
                    br_vec_.append((br_vec[i], data[br_vec[i]]))
         
            # далее идет поиск "заветной" скобки
            # поиск осуществляется по средством перебора массива со скобками, игнорируя найденную
            # первую скобку для того, чтобы определить последнюю скобку используется коэффициент, 
            # который увеличивается на один, если встречается открывающая (левая) скобка и 
            # уменьшается на единицу, если встречается закрывающая скобка
            # в конце каждой итерации (проверки очередной скобки) проверяется на найдена ли нужная скобка
            # с помощью коэффициента, который в таком случае будет равняться нулю
            br_it = 1
            end_br_pos = None
            for i in br_vec_:
                # так как алгоритм выше ищет всех скобки в строке, необходимо сразу отсеять скобки,
                # находящиеся до первой скобки (связанной с рассматриваемым корнем)
                if i[0] <= begin_br_pos:
                    continue
                if (i[1] == '('):
                    br_it += 1
                else:
                    br_it -= 1
                if (br_it == 0):
                    end_br_pos = i[0]
                    break
                    
            # так как координаты скобок корня найдены, то можно заменить корень на степень
            data = data[:end_br_pos + 1] + "**0.5 " + data[end_br_pos + 1:]
            
            # удаление слова 'sqrt'
            data = data[:pos] + data[pos+4:]
   
        # если входное выражение (или уже преобразованное) не содержит корней
        # возвращаем его
        else:
            return data

def swap_sqrt_fast(data:str):
    '''
        Функция, копирующая первоначальную swap_sqrt(data). Является немного эффективней, так 
        отсутствуют ненужные шаги.
        
        Функция, получающая на вход строку (любое математическое выражение) преобразующее 
        все корни (sqrt(...)) в степенной вид (**0.5). 
        
        Пример. 
        Вызов функции с следующим аргументом: swap_sqrt('sqrt(3) / sqrt(4)').
        Вернет: '(3)**0.5 / (4)**0.5'
        
        Сложность данной задачи заключается в корректном трекинге скобок. Так как существуют
        усложненные выражения типа: 'sqrt(sqrt(sqrt())) - sqrt()' и похожее. Так как скробок 
        может быть очень много, необходимо неати именно ту, что отвечает за корень и добавить
        после нее степень.
        
        Найденное решение не является эффективным.
    '''
    
    # цикл while необходим для преобразования выражений, содежащих несколько корней
    # для каждого 'sqrt' будет проходит одна итерация цикла while
    while (True):
        # проверяем, содержится ли в входном выражении корни
        if ('sqrt' in data):
            
            # находим первый корень (поиск ведется до первого встреченного значения слева строки)
            pos = data.find('sqrt')
            
            # координата левой скобки, связанной с корнем
            # литера 4 - есть длина слово sqrt
            begin_br_pos = pos + 4

            # объявляем два вектора, которые будут содержать координаты (позиции)
            # левых и правых скобок в входной строке, что важно - до изменения
            left_br_poses = []
            right_br_poses = []
    
            # вектор, который будет хранить позиции и варианты скобок
            br_vec_ = []

            # просматриваем каждый символ строки и записываем в вектор позицию и значение скобки
            for i in range(0, len(data)):
                if (data[i] == '('):
                    br_vec_.append((i, '('))
                if (data[i] == ')'):
                    br_vec_.append((i, ')'))

            # далее идет поиск "заветной" скобки
            # поиск осуществляется по средством перебора массива со скобками, игнорируя найденную
            # первую скобку для того, чтобы определить последнюю скобку используется коэффициент, 
            # который увеличивается на один, если встречается открывающая (левая) скобка и 
            # уменьшается на единицу, если встречается закрывающая скобка
            # в конце каждой итерации (проверки очередной скобки) проверяется на найдена ли нужная скобка
            # с помощью коэффициента, который в таком случае будет равняться нулю
            br_it = 1
            end_br_pos = None
            for i in br_vec_:
                # так как алгоритм выше ищет всех скобки в строке, необходимо сразу отсеять скобки,
                # находящиеся до первой скобки (связанной с рассматриваемым корнем)
                if i[0] <= begin_br_pos:
                    continue
                if (i[1] == '('):
                    br_it += 1
                else:
                    br_it -= 1
                if (br_it == 0):
                    end_br_pos = i[0]
                    break
                    
            # так как координаты скобок корня найдены, то можно заменить корень на степень
            data = data[:end_br_pos + 1] + "**0.5 " + data[end_br_pos + 1:]
            
            # удаление слова 'sqrt'
            data = data[:pos] + data[pos+4:]
   
        # если входное выражение (или уже преобразованное) не содержит корней
        # возвращаем его
        else:
            return data


def solve_(leftpart, rightpart, variable):
    lex = sp.parse_expr(leftpart)
    rex = sp.parse_expr(rightpart)
    
    eq = sm.Eq(lex, rex)
    
    ans, = sm.solve(eq, 'x')
    
    return(ans)

def solve__(exp, parameters, variable, list_mode=False):
    '''
        Функция созданная с целью соединить решение уравнений из Sympy совместно с единицами измерения
        из pint.
        Минусом является - небольшая гибкость входные данных.
        
        Функция на вход принимает 3 обязательных и 1 необязательный (флаг) аргументов.
        
        Первый входной параметр есть список уравнений (возможно одно уравнение), записанный в виде tuple
        ("a * x - c * y",   '22') - где левая часть является левой частью уравнения, а правая часть - правой. лол
        Пример:
        exp = [
            ("a * x - c * y",   '22'),
            ('b*x + d',         '45')
        ]
        
        Второй входной параметр есть список известных параметров. Конечно, можно заранее заменить 
        переменные значениями, но функция способна сделать это by itself.
        Значения подаются так же в виде списка, состоящего из tuple's, левая часть которых есть
        строчное обозначение переменной в уравнениях из exp, правая же - их значение, в виде int, float, pint.Quantity
        Пример:
        parameters = [
            ('a', aVALUE),
            ('b', aVALUE),
            ('c', aVALUE),
            ('d', aVALUE),
        ]
        
        Третий параметр есть список необходимых для нахождения величин. Так же записываются в виде списка, но без tuple, 
        в строчном виде (так же как и записаны в exp).
        Пример:
        variable = [        variable = [
            'x',                ('x'),
            'y'                 ('y'),
        ]                   ]
        
        Четвертный параметр есть флаг, который устанавливает ражим вывода ответа. Данный флаг имеет смысл использовать, 
        если уравнения могут дать не одно решение, а несколько, как если бы мы находили решение квадратного уравнения.
        Если list_mode = False - решения будут выводиться вложенными, то есть один вариант решения (найденные значения
        переменных, описанных в variable) в собственном списке, вложенном в общий список.
        пример вывода:
        solve__(...) = [
            [10 meter, 30 kg/m],
            [-10 meter, 40 kg/m]
        ]
        Если list_mode = True - решения вернуться в одном списке:
        solve__(...) = [
            10 meter, 
            30 kg/m,
            -10 meter, 
            40 kg/m
        ]
        
    '''
    # чтобы было. 
    if ((isinstance(exp, list))and(isinstance(variable, list))):
        # если неизвестных столько же сколько и уравнений
        if (len(exp) == len(variable)):
            '''
                exp = [
                    ("a * x - c * y",   '22'),
                    ('b*x + d',         '45')
                ]
                parameters = [
                    ('a', aVALUE),
                    ('b', aVALUE),
                    ('c', aVALUE),
                    ('d', aVALUE),
                ]
                variable = [
                    ('x'),
                    ('y')
                ]
            '''                
            # 1. формируем вектор уравнений
            # то что будет хранить уравнения
            eq_vec = []
            # для каждой строчки в входном вектора
            for it in exp:
                # левая часть итерируемой строчки
                lval = sp.parse_expr(it[0])
                # правая часть итерируемой строчки
                rval = sp.parse_expr(it[1])
                # добавляем в новому вектору, перенося на одну сторону
                # чтобы достичь вида: f(...) = 0
                eq_vec.append(lval - rval)
            
            # 2. формируем ветор решений
            sol = sm.solve(eq_vec, variable, dict=True)
            
            #print('sol:',sol)
            
            # 3. формируем численное решение
            param_vec = []
            # преобразуем значения из pint в вид, необходимый для sympy
            for it in parameters: 
                param_vec.append((it[0], make_correct_str_base(it[1])))

            # вектор, который будет хранить выходные данные (решение системы или уравнения)
            output_values = []
            
            # если необходимо вывести все в одном списке
            if (list_mode):
                # чтобы было
                if (isinstance(sol, list)):
                    # для каждого решения системы:
                    for j in sol:
                        # для каждой находимой переменной необходимо заменить общий вид на числовые значения
                        for i in range(0, len(variable)):
                            '''
                                Замена переменных на значения происходит в несколько этапов.
                                1. Один известный параметр (входной) в одну итерацию i. Так
                                сначала получаем название этой переменной.
                                2. с помощью sm.symbols преобразовываем переменную в символьную
                                переменную SymPy, чтобы можно было ее отыскать по словарю, который
                                получается при использовании функции SymPy.solve (определяемые из системы 
                                переменные: значения в символном виде)
                                3. получаения с помощью get переменной, которую будем заменять в ланной итерации i
                                4. далее с помощью функции SymPy.subs заменяем одну значение переменной на ее значение
                                
                                Так, проделывая все это с каждой переменной, получаем числовое значение.
                            '''
                            str_version = str(j.get(sm.symbols(variable[i])).subs(param_vec))
                            str_version = str_version.replace('-', '+(-1)*')

                        
                            # если были полученые корни в решении - их необходимо заменить на степени, так как
                            # pint не работает с корнями
                            if ('sqrt' in str_version):
                                str_version = swap_sqrt(str_version)

                            # Преобразуем в pint.Qunatity и добавляем к выходному вектору значений
                            output_values.append(u(str_version))

                    return (output_values)
            
            # если необходимо вывести вложенно
            else:
                if (isinstance(sol, list)):
                    for j in sol:
                        output_pack = []
                        for i in range(0, len(variable)):
                            str_version = str(j.get(sm.symbols(variable[i])).subs(param_vec))
                            str_version = str_version.replace('-', '+(-1)*')
                            #print("str_verstion:",str_version)
                        
                            if ('sqrt' in str_version):
                                str_version = swap_sqrt(str_version)
                        
                            output_pack.append(u(str_version))
                        output_values.append(output_pack)
                    return (output_values)
                   


                 
  

# sol = (solve__(
#     [
#        ('4 * l_g**2 - 2 * y_d * v_2t/v_2z', '22 * meter**3 / kg'),
#        ('6*x + d_k', '45 * meter**3 / kg')
#     ],
    
#     [
#         ('y_d', 10 * u('meter**3 / kg')),
#         ('x', 33 * u('meter**3 / kg')),
#         ('v_2t', 220 * u('meter**3 / kg')),
#         ('v_2z', 330 * u('meter**3/ kg'))
#     ],
    
#     [
#         ('l_g'),('d_k')
#     ],
#     list_mode=False
# ))
# for i in sol:
#     print(i)

#lstr = "sqrt(d*(10+4)/(34 * k))" 
#llstr = "sqrt(sqrt(a/b*(a+b)))"
#lllstr = '(sqrt((sqrt((a**2+b)*100*c-1488*d+3228))/(0.5**0.5 + 0.25*c - 3826*d))/(sqrt(35.5*a*b)*sqrt(3561*d+a*b))'
#llllstr = 'sqrt(44)/sqrt(55)(10 + b)/(sqrt(sqrt(33+90)*10))'
#lllllstr = 'sqrt(sqrt(44) * 10 )'
#print(llstr)
# print(swap_sqrt_fast(llllstr))

#abc = [10, 33, 2, 1]
#print(abc)
#abc.sort()
#print(abc)
#print("sort", abc)


def sub_(exp, data):
    #if (isinstance(data, dict)):    # {}
    
    str_data = [] 
    if (isinstance(data, list)):    # []
        
        for it in data:
            str_data.append((it[0], make_correct_str_base(it[1])))
        for it in data:
            exp = exp.subs(str_data)
            
    return u(str(exp))

def return_positive(first, second):
    if (first > 0):
        return first
    else:
        return second

# print(solve__(
#     [
#         ('l_2z**2 + l_2z*d_k','l_21*d_1*(v_2z/v_2t)'),
#     ],
    
#     [
#         ('d_k', pint.Quantity(0.7838715992465253, "meter")),
#         ('l_21', pint.Quantity(0.01612840075347477, "meter")),
#         ('d_1', pint.Quantity(0.9, "meter")),
#         ('v_2z', pint.Quantity(0.155, 'meter**3/kg')),
#         ('v_2t', pint.Quantity(0.042, 'meter**3/kg'))
#     ],
    
#     [
#         'l_2z'
#     ],
#     list_mode=True
# ))
# '-0.4551735955658758*kilogram*(meter**8/kilogram**2)**0.5 /meter**3 - 0.3919357996232627*meter'
#print(u('+(-1)0.4551735955658758*kilogram*(meter**8/kilogram**2)**0.5 / (meter**3) +(-1)* 0.3919357996232627*meter'))
#print(u('-0.4551735955658758*kilogram*(meter**8/kilogram**2)**0.5 / (meter**3) -0.3919357996232627*meter'))
# exp33 = [
#     ('d_k', d_k),
#     ('l_21', l_21),
#     ('v_2z', v_2z)
# ]
# print(sub_(solve_('d_k / l_21 * x','v_2z', 'x'), exp33))
#     if (isinstance(data, tuple)):   # ()
            

# exp_1l = "4 * x**2 - 2 * y"
# exp_1r = '22'
# exp_2l = '6*x + y'
# exp_2r = '45'

# #print(sp.parse_expr(str("4 * x - 2 * y")))
# exp_1l_ = sp.parse_expr(exp_1l)
# exp_1r_ = sp.parse_expr(exp_1r)
# exp_2l_ = sp.parse_expr(exp_2l)
# exp_2r_ = sp.parse_expr(exp_2r)

# eq_1 = exp_1l_ - exp_1r_#sm.Eq(exp_1l_, exp_1r_)
# eq_2 = exp_2l_ - exp_2r_#sm.Eq(exp_2l_, exp_2r_)
# #print(eq_1)
# ans_ = sm.solve(['4 * x - 2 * y - 22', '6*x + y - 45'], ['x','y'], dict=True)
# print(ans_)
# ans_ = sm.solve(['4 * x**2 - 2 * y**2 - 22', '6*x + y - 45'], ['x','y'], dict=True)
# print(ans_)
#print(ans_)
#print(solve_('d_k / l_21 * x','v_2z', 'x'))
#ans.subs([('l_21', make_correct_str_base(l_21)), ('v_2z', make_correct_str_base(l_21)), ('d_k', make_correct_str_base(l_21))])
#ans_next = ans.subs([('l_21', make_correct_str_base(l_21)), ('v_2z', make_correct_str_base(l_21)), ('d_k', make_correct_str_base(l_21))])
#ans_pint = u(str(ans_next))
#print(ans_pint)
#print(type(ans_pint))

#def names(data):
    #print(varname.nameof(data))
    
#names(d_k)

#print(varname.nameof(some_variable))

#eq2 = sm.Eq(d_k/l_21*x,v_2t)
#print(d_k/l_21*x)

#eq1 = sm.Eq(sp.parse_expr(make_correct_str_base(s1) + " + " + make_correct_str_base(s2), sp.parse_expr('x')))

#print(sp.parse_expr(exp1))
#print(sm.solveset(eq1, sp.parse_expr('x')))
#print(sm.Eq(sp.parse_expr(make_correct_str_base(s1)), sp.parse_expr('x')))
#sm.nsolve(sm.Eq(sp.parse_expr(make_correct_str_base(s1),sp.parse_expr('x')),sp.parse_expr('x')))


#eq1 = sm.Eq((l_2z)**2  + l_2z * d_k, l_21*d_21*v_2z/v_2t)
#print(eq1)
#ans = sm.nsolve(eq1, l_2z)
#print(ans)

