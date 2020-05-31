function alpha = armijo(x,objFunc,q)


% 1 obtain search direction
d = -q;
obj = objFunc(x);
sigma = .01;
beta = .5;              
alpha = 1;
newobj = objFunc(x + alpha*d);
nf = 1;

while (newobj-obj)/alpha > sigma*q'*d
  alpha = alpha*beta;
  newobj = objFunc(x+ alpha*d);
  nf = nf+1;
end


