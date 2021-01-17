import sys

sys.path.append("..")



# import matplotlib.pyplot as plt


from sympy import *

# import numpy as np


# from colors.mycolors import *


init_printing(use_unicode=True, wrap_line=False)




x     = Symbol('x',  real=True)

mu    = Symbol('mu',    real=True)
alpha = Symbol('alpha', real=True, positive=True)
nu    = Symbol('nu',    real=True, positive=True)
beta  = Symbol('beta',  real=True, positive=True)


Delta = 3 * sqrt(2) / (alpha * nu)

theta = (1/4) * log((1+nu)/(1-nu))



u = (alpha*nu/3) * (tanh(x/Delta+theta) - tanh(x/Delta-theta))

u    = simplify(u)

# print(u)

# u_xx = diff(u, x, x)

u_x = Derivative(u, x)   

# u_x = simplify(u_x)

print(u_x)




