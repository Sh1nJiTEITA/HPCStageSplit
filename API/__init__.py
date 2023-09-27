
from .snj_solvers import *
from .TCPv2 import *

from .HPC_SPLIT import *
from .LATEX_OUTPUT import *

from .losses import *
from .HPC_SPLIT import *

from .CustomBladesProfiles import *
from .Stage import *

from .lstr import *

from .tlstr import tlstr
from .losses import calculate_losses


#make_math_string("3123", a=10, b=30)

GLOBAL_DF = pd.DataFrame(columns=[
        'n', 
        'rho_k', 
                                   
        'd_1_1', 
        'd_next_1', 
        'alpha_1eef_1',
        'Z_1',
        'Z_new_1',
        'Z_ratio_1',
        'q_t_1',
        'Delta_H_1',
        'H_ave_1',
        'l_11_1',
        'd_k_1',
        'theta_1',
        
        
        'd_1_2', 
        'd_next_2', 
        'alpha_1eef_2',
        'Z_2',
        'Z_new_2',
        'Z_ratio_2',
        'q_t_2',
        'Delta_H_2',
        'H_ave_2',
        'l_11_2',
        'd_k_2',
        'theta_2',
        
        'd_1_3',
        'd_next_3', 
        'alpha_1eef_3',
        'Z_3',
        'Z_new_3',
        'Z_ratio_3',
        'q_t_3',
        'Delta_H_3',
        'H_ave_3',
        'l_11_3',
        'd_k_3',
        'theta_3'
        ])

print("API loaded")