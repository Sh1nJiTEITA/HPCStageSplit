# import sympy
# from pint import get_application_registry
# from pint import Quantity
# from sympy import sympify
# Quantity._sympy_ = lambda s: sympify(f'{s.m}*{s.u:~}')
# un = get_application_registry()

# x = sympy.symbols('x')

# y = 10 * un('meter^2')

# eq1 = sympy.Eq(x*4 + y, 0)

# f, = sympy.solveset(eq1, x)

# #print(f)

# #print(type(f))

# print(un(str(f)))

# #print(sympy.solveset(eq1, x))

# Import NumPy
import numpy as np

# Disable Pint's old fallback behavior (must come before importing Pint)
# import os
# os.environ['PINT_ARRAY_PROTOCOL_FALLBACK'] = "0"

# Import Pint
import pint
ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

# # Silence NEP 18 warning
# import warnings
# with warnings.catch_warnings():
#     warnings.simplefilter("ignore")
#     Q_([])
    
legs1 = Q_(np.array([10, 20, 30, 40]), 'meter')
print(legs1)