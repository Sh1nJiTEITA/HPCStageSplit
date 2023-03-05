
import sys, os

import inspect 
loc_path = inspect.stack()[0][1]

loc_path = loc_path[:loc_path.rfind("\\")]
#print(loc_path)
loc_path = loc_path[:loc_path.rfind("\\")]
#print(loc_path)
loc_path = loc_path[:loc_path.rfind("\\")]
#print(loc_path)

sys.path.append(loc_path + '\\API')

# for i in sys.path:
#     print(i)
    
import TCPv2
import HPC_SPLIT as spl
import snj_solvers
from LATEX_OUTPUT import display_latex_ex_2 as _dl2

print("Files loaded")
