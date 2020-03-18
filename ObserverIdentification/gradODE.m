function grad = gradODE(t,x,X1,nrow, ncol,C)


U = reshape(x(1:nrow*ncol),nrow,ncol);
V = U*U'; n = nrow/2;
P1 = V(1:n,1:n) ;
Q1 = V(1:n,n+1:end);
P2 = V(n+1:end,n+1:end) ;
Q2 = V(n+1:end,1:n);
L = x(nrow*ncol + 1:end);
acc = 100;


grad1 = -acc*[(P1*X1 - Q1 + L*C*P1)*(X1 + L*C)', (-X1*P2 + Q1 - L*C*P2) ;...
            (-X1*P1 + Q2' -L*C*P1)'  , (P2*X1 - Q2' + L*C*P2)*(X1 +L*C)']*U ;
grad2 = -acc*(P1*X1 - Q1 + L*C*P1)*(C*P1)';

grad11 = reshape(grad1,nrow*ncol,1);
grad22 = reshape(grad2,n,1);

grad = [grad11;grad22];


end 