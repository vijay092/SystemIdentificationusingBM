function f = CostFn(x,X1,X2, X3, X4, nrow,ncol);

U = reshape(x, nrow, ncol);
V = U*U'; 
n = nrow/2;
f = (norm(X1*V(1:n,1:n) - V(1:n,n+1:end),'fro'));

end









