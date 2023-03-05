from IPython.display import Latex, display
import numpy as np

def output_value_var_str(var, value):
    if not((isinstance(value, int))or
           (isinstance(value, float))or
           (isinstance(value,np.float64))):
        value = "{}\ {}".format(value.magnitude, value.units)
    #return (f"""\\begin{{align}}{var} = {value}\\end{{align}}""")
    return (f"""{var} = {value}""")
    #  return Latex(f"""\\begin{{align}}
                 #\\{var} = {value}\\end{{align}}""")
def op_latex_output(data):
    begin_str = f'\\begin{{dcases}}'
    end_str = f'\\end{{dcases}}'
    output_str = ""
        
    for i in data:
        
        output_str+=(output_value_var_str(i[0], i[1]))
        output_str+=r" \\ "

    re_str = begin_str + output_str
    re_str += end_str
    
    return re_str


def display_latex_ex(data):
    display(Latex(op_latex_output(data)))

