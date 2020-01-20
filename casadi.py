

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:02:47 2020

@author: sanja
"""
from casadi import *


x = SX.sym('x'); 
# nlp = {'x':x, 'f':x.T @ x, 'g':0}
# S = nlpsol('S', 'ipopt', nlp)
# print(S)


# x0 = [1,1]

# r = np.zeros(2)
# r = S(x0,\
#       lbg=r, ubg=r)
# x_opt = r['x']
# print('x_opt: ', x_opt)

