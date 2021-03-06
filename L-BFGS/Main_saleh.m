% ======================================================
% This problem identifies an observer directly from 
% data. We solve the optimization problem using 
% L-BFGS


clc 
clear all

% Generate a random state-space system
n = 2;
csys = rss(n); 
ts = 0.1;
sys = c2d(csys, ts);


% Generate time-series data
t = 0:ts:3;
X = zeros(n,length(t));
Y = zeros(1,length(t)-1);
A = sys.A; B = sys.B; C = sys.C; D = sys.D;
n = length(A);
for i = 1:(length(t)-1)
    X(:,i+1) = A*X(:,i) + B;
    Y(:,i) = C*X(:,i) + D;  
end 
figure(1);
plot(Y,'Linewidth',2);

% Construct snapshots:
X1_til = X(:,1:end-1);
X2_til = X(:,2:end);
u = ones(length(t),1)';
U1 = u(1:end-1);
X = X2_til*Rinv([X1_til;U1]);
X1 = X(1:n,1:n);

% Sized of U
nrow = 2*n;
ncol = 2*n;


% specify input
mem          = 10;        % number of past gradients and function values used for inverse hessian contruction
x            = rand(nrow*ncol + n,mem);
x(:,1)       = rand(nrow*ncol + n,1);  % starting value for parameters
objFuncValue = rand(1,mem);
dx           = rand(nrow*ncol + n,mem);

s_k = ones(nrow*ncol + n,mem-1);
y_k = ones(nrow*ncol + n,mem-1);
r_k = ones(mem-1,1);
a_k = ones(1,mem-1);
       
% 1st calculation of objective function and gradient

objFunc         = @(x) ObjFun(x,X1,nrow, ncol,C);
objFuncValue(1) = objFunc(x(:,1));
objFuncValue(2) = objFuncValue(1) + 1;
dx(:,1)          = numDiff(objFunc,x(:,1));

% iterate
iter      = 0;
numOfIter = 1000 ;
prec      = 1e-9;
tic

% convergence if gradient smaller than prec, change in objective function
% smaller than prec or maximum number of iteration reached...
while iter < numOfIter %&& abs((objFuncValue(2)-objFuncValue(1))/objFuncValue(1))>prec && norm(dx(:,1))>prec
    
     % inverse hessian update
    q = dx(:,1);
    for i = 1:min(iter,mem-1)
        a_k(i) = r_k(i)*s_k(:,i)'*q;
        q      = q - a_k(i)*y_k(:,i);
    end
    z = s_k(:,1)'*y_k(:,1)/(y_k(:,1)'*y_k(:,1))*q; %approxm of (H*dx)
    for i = min(iter,mem-1):-1:1
        b = r_k(i)*y_k(:,i)'*z;
        z = z + s_k(:,i)*(a_k(i)-b);
    end
   
    
    % 1 obtain search direction
    dir = 1e4*-z;
    
    % 2.1 linesearch to to find acceptable stepsize alpha
    % you can replace this part with the method of choice for caluclating
    % step size, Armijo, backtracking, etc. For CFD calculations this is an
    % open question (in terms of technicality; the fudamentals of math is
    % known).
    % should change it for consistency in size; let's stick to backtracking
    % for now
    %alpha = nocLineSearch(objFunc,@(x) numDiff(objFunc,x),x,dir,dir'*dx,objFuncValue);
    alpha = backtrackingLineSearch(objFunc,objFuncValue(1),x(:,1),dx(:,1),dir);
    %alpha = alpha/2;

    % 3
    p = alpha*dir;

        
    % 2.2 update x
    x(:,2:end) = x(:,1:end-1);
    x(:,1)     = x(:,1) + p;
    s_k        = -diff(x,[],2);
    
    % update objective function (and remember old objective function to
    % check for convergence)
    objFuncValue(2) = objFuncValue(1);
    objFuncValue(1) = objFunc(x(:,1));
    
    % increment iteration counter
    iter = iter + 1;
    
    fprintf(1,'Iteration %d: alpha=%f, OF=%f\n',iter,alpha,objFuncValue(1));
    
    % 4 calculate difference of gradients
    dx(:,2:end) = dx(:,1:end-1);
    dx(:,1)     = numDiff(objFunc,x(:,1));
    y_k         = -diff(dx,[],2);

    r_k = 1./diag(y_k'*s_k);
    
    
    
end
toc

lenLMI = nrow*ncol;
lenL = nrow/2;
% Get back required matrices:
U_opt = reshape(x(1:lenLMI,1),nrow, ncol) ;
L_opt = x(lenLMI+1:end,1);
V_opt = U_opt * U_opt';
P_opt = V_opt(1:n,1:n);
Q_opt = V_opt(1:n,n+1:end) + L_opt*C*P_opt;
A_opt = (Q_opt * inv(P_opt));
B_opt = X(:,n+1:end); 


% The output equation:
Rx = Rinv([X1_til;U1])*blkdiag(P_opt, 1); 
Y = Y*Rx*Rinv([X1_til;U1]*Rx);
C_opt = Y(1:n);
D_opt = Y(n+1:end);

x_id = zeros(n,length(t));
y_id = zeros(1,length(t)-1);

for i = 1:(length(t)-1)
    x_id(:,i+1) = A_opt*x_id(:,i) + B_opt;
    y_id(:,i) = C_opt*x_id(:,i) + D_opt;
end
figure(1);
hold on;
plot(y_id,'r*','Linewidth',2)