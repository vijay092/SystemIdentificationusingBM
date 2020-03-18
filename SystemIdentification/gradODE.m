function grad = gradODE(t,x,X1,nrow, ncol)

% Your x is your U.
% it is going to be a vector.


% shape


U = reshape(x,nrow,ncol);
V = U*U'; n = nrow/2;
P1 = V(1:n,1:n) ;
Q1 = V(1:n,n+1:end);
P2 = V(n+1:end,n+1:end) ;
Q2 = V(n+1:end,1:n);


grad1 = -100*[(X1'*X1*P1 - Q1'*X1), (-X1*P2 + Q1) ;...
            (-X1*P1 + Q2')'  , (X1'*X1*P2 - Q2*X1)]*U ;

grad = reshape(grad1,nrow*ncol,1);



end 