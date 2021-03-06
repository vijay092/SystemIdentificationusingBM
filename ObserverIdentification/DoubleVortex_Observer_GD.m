% =========================================================================
%% This code identifies an observer, L purely using data. 
% The SDP formulation
% of the problem described in overleaf (Observer Identification using BM)
% is not scalable. This code uses
% Burer-Monteiro to eliminate the LMI, translate the problem into an
% unconstrained non-convex problem. 

% =========================================================================
% Author: Sanjana Vijayshankar
% =========================================================================
clear all; close all

% Generate a random state-space system
n = 5;
csys = rss(n); 
ts = .1;
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
lenLMI = nrow*ncol;
lenL = nrow/2;

%% Gradient descent with armijo step size rule
iterMax = 100;
x = rand(lenLMI + lenL ,1);
gradFunc = @(x) gradODE(t,x,X1,nrow, ncol,C)
objFunc = @(x) ObjFun(x,X1,nrow, ncol,C)
iter = 0
while iter < iterMax
    grad = gradFunc(x);
    alpha = armijo(x,objFunc,grad);
    x = x - alpha*grad;
    iter = iter + 1;
    fprintf('Iteration = %d, \t alpha = %f,\t Obj Value = %f \n',iter,alpha,objFunc(x))
end

x_opt = x;
% Get back required matrices:
U_opt = reshape(x_opt(1:lenLMI),nrow, ncol) ;
L_opt = x_opt(lenLMI+1:end);
V_opt = U_opt * U_opt';
P_opt = V_opt(1:n,1:n);
Q_opt = V_opt(1:n,n+1:end) - L_opt*C*P_opt;
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

% To make sure the observer will converge..
fprintf('Maximum eigen value of A+LC = %f \n',max(eig(A_opt + L_opt*C_opt)));
