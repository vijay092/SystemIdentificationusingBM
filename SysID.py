# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 14:03:00 2020

@author: sanja
"""


'''

This code simulates a system identification problem cast as an LMI. 
The idea is driven by a paper by Lacy and Berstein

"Subspace identification with guaranteed stability using constrained optimization"

This paper discusses a stability preserving framework for system identification 
 of linear systems using LMIs 

Going to implement this using picos 

'''


import picos as pic
import cvxopt as cvx
import numpy as np



A = np.array([[1,0],[4,5]])

A_cvx = cvx.matrix(A)

