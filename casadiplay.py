

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:02:47 2020

@author: sanja
"""
from casadi import *


x = SX.sym('x',2,2); 
z = x @ x.T
f = reshape(x,4,1)
nlp = {'x':f, 'f':trace(z), 'g':x[0]}
S = nlpsol('S', 'ipopt', nlp)
print(S)


r = S(x0=[100,1,3,4],lbg=[1], ubg=[1])
x_opt = r['x']
print('x_opt: ', x_opt)


import cvxpy as cp
import numpy

# Problem data.
m = 30
n = 20
numpy.random.seed(1)
A = numpy.random.randn(m, n)
b = numpy.random.randn(m)

# Construct the problem.
x = cp.Variable(n)
objective = cp.Minimize(cp.sum_squares(A*x - b))
constraints = [0 <= x, x <= 1]
prob = cp.Problem(objective, constraints)

# The optimal objective is returned by prob.solve().
result = prob.solve()
# The optimal value for x is stored in x.value.
print(x.value)
# The optimal Lagrange multiplier for a constraint
# is stored in constraint.dual_value.
print(constraints[0].dual_value)