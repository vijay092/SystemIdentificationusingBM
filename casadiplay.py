

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:02:47 2020

@author: sanja
"""
from casadi import *


x = SX.sym('x',2,2); 
A = 3*np.eye(2);
z = A @ x @ x.T
f = reshape(x,4,1)
nlp = {'x':f, 'f':trace(z) , 'g':trace(A @ x)}
S = nlpsol('S', 'ipopt', nlp)
print(S)


r = S(x0=[100,1,93,4],lbg=[1], ubg=[1])
x_opt = r['x']
print('x_opt: ', x_opt)

