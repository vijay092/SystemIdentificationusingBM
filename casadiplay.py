

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:02:47 2020

@author: sanja
"""
from casadi import *


x = SX.sym('x',2); 
b = [1]
nlp = {'x':x, 'f':x.T @ x + b, 'g':x}
S = nlpsol('S', 'ipopt', nlp)
print(S)


r = S(x0=[0,0],lbg=[1,1], ubg=[1,1])
x_opt = r['x']
print('x_opt: ', x_opt)