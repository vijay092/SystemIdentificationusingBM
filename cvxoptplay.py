# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 12:04:15 2020

@author: sanja
"""


from cvxopt import matrix, solvers

c = matrix([10., 1.,10.])
G =[ matrix([[1., 0., 0., 0.], [0., 1., 1., 0.],[0., 0., 0., 1.]])]
h = [ matrix([[0., 0.], [0., 0.]]) ]
sol = solvers.sdp(c, Gs=G, hs=h)

print(sol['x'])


    