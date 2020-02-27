function grad = gradODE(t,x,X1)

% Your x is your U.
% it is going to be a vector.


% shape
nu = sqrt(length(x)); n = nu/2;

U = reshape(x,nu,nu);
V = U*U';
P = V(1:n,1:n);
Q = V(1:n,n+1:end);
t
grad = -100*[(X1'*X1*P - Q'*X1), (-X1*P + Q) ;...
            (-X1*P + Q)'  , (X1'*X1*P - Q'*X1)]*U; 

grad = reshape(grad,nu^2,1);



end 