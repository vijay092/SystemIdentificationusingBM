# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 12:04:15 2020

@author: sanja
"""


from cvxopt import matrix, solvers
import picos as pic 
import cvxopt as cvx


'''
I'm trying to munderstand this cvxopt. I want to be able to solve an sdp using
cvxopt. 


'''


c = matrix([1., 1.,1.])
G =[ matrix([[-1., 0., 0., 0.], [0., -1., -1., 0.],[0., 0., 0., -1.]])]
h = [ matrix([[1., 1.], [1., 0.]]) ]
sol = solvers.sdp(c, Gs=G, hs=h)

print(sol['x'])

#---------------------------------#
# First generate some data :      #
#       _ a list of 8 matrices A  #
#       _ a vector c              #
#---------------------------------#
A=[ cvx.matrix([[1,0,0,0,0],
                [0,3,0,0,0],
                [0,0,1,0,0]]),
cvx.matrix([[0,0,2,0,0],
                [0,1,0,0,0],
                [0,0,0,1,0]]),
cvx.matrix([[0,0,0,2,0],
                [4,0,0,0,0],
                [0,0,1,0,0]]),
cvx.matrix([[1,0,0,0,0],
                [0,0,2,0,0],
                [0,0,0,0,4]]),
cvx.matrix([[1,0,2,0,0],
                [0,3,0,1,2],
                [0,0,1,2,0]]),
cvx.matrix([[0,1,1,1,0],
                [0,3,0,1,0],
                [0,0,2,2,0]]),
cvx.matrix([[1,2,0,0,0],
                [0,3,3,0,5],
                [1,0,0,2,0]]),
cvx.matrix([[1,0,3,0,1],
                [0,3,2,0,0],
                [1,0,0,2,0]])
]

c = cvx.matrix([1,2,3,4,5])

#create the problem, variables and params
prob_SDP_c_primal=pic.Problem()
AA=[cvx.sparse(a,tc='d') for a in A] #each AA[i].T is a 3 x 5 observation matrix
s=len(AA)
AA=pic.new_param('A',AA)
cc=pic.new_param('c',c)
mu=prob_SDP_c_primal.add_variable('mu',s)

#define the constraints and objective function
prob_SDP_c_primal.add_constraint(
        pic.sum(
        [mu[i]*AA[i]*AA[i].T for i in range(s)], #summands
        'i', #index
        '[s]' #set to which the index belongs
        )
        >> cc*cc.T )
prob_SDP_c_primal.add_constraint(mu>0)
prob_SDP_c_primal.set_objective('min',1|mu)

#solve the problem and retrieve the weights of the optimal design
print(prob_SDP_c_primal)
prob_SDP_c_primal.solve(verbose=0,solver='cvxopt')
w=mu.value
w=w/sum(w) #normalize mu to get the optimal weights
print()
print('The optimal design is:')
print(w)