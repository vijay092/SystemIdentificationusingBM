

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:02:47 2020

@author: sanja
"""
from casadi import *


x = SX.sym('x',2,2); 
z = x.T @ x
f = reshape(x,4,1)
nlp = {'x':f, 'f':z[0] + z[3], 'g':z[0]}
S = nlpsol('S', 'ipopt', nlp)
print(S)


r = S(x0=[0,1,3,4],lbg=[1], ubg=[1])
x_opt = r['x']
print('x_opt: ', x_opt)