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
from scipy import signal
import numpy.linalg as la

# Function to compute R inverse
def Rinv(v):
    vt = np.transpose(v)
    r = np.dot(vt,la.pinv(np.dot(v,vt)));
    return r


# State Matrices of the system to be identified
A = np.array([[-1, 0],[0, -2]])
B = np.array([[0],[1]])
C = np.array([0,1])
n = len(A)

# Define the system
system = signal.lti(A, B, C, 0.)
t = np.linspace(0, 5)   
u = np.ones_like(t)
tout, y, x = signal.lsim(system, u, t)
x = np.transpose(x)

# Construct relevant matrices 
X1_til = x[:-1];
X2_til = x[1:];
U1 = u[:-1]
XU =  Rinv(np.vstack((X1_til,U1)));
X = np.dot(X2_til, XU)
X1 = X[:n,:n]
X2 = X[:,n:]


# Parameters for the SDP
dell = 1e-2;
eps = 25;

# The problem in picos 
sysid = pic.Problem();
X1_cvx = cvx.matrix(X1);
X1_para=pic.new_param('X1p',X1_cvx)

# The arguments
P = sysid.add_variable('P',(2,2));
Q = sysid.add_variable('Q',(2,2));

# Constraint
sysid.remove_all_constraints()
sysid.add_constraint((P & Q.T) // (Q & P) >> 0 )

# Objective fucnction

