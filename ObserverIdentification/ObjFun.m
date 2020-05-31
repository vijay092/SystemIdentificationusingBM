function F = ObjFun(x,X1,nrow, ncol,C)


U = reshape(x(1:nrow*ncol),nrow,ncol);
V = U*U'; n = nrow/2;
P1 = V(1:n,1:n) ;
Q1 = V(1:n,n+1:end);
P2 = V(n+1:end,n+1:end) ;
Q2 = V(n+1:end,1:n);
L = x(nrow*ncol + 1:end);

F = norm(X1*P1 - Q1 + L*C*P1,'fro')^2;% + norm(X1*P2 - Q2 - L*C*P2,'fro')^2;



end