% ======================================================
% This problem identifies an observer directly from 
% data. We solve the optimization problem using 
% L-BFGS


clear all
close all
load StateMatricesDoubleVortex
load x_opt

%csys = rss(n); 
ts = 2;
idx = 1:1000;
% select the number of sensors:

%sys = c2d(csys, ts);
A = A(idx,idx);
b = b(idx);
n = length(A);
ny = 1;
% Generate time-series data
t = 0:ts:70;
x0 = zeros(n,1) ;
[t,X_ts] = ode45(@(t,x) A*x + b, t, x0);
X = X_ts';
y_idx = 1:ny;

C = zeros(ny,n);
sens = 1;
C(sens) = 1;

Y = C*X(:,1:end-1);

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

options = struct('GradObj','on','Display',...
'iter','LargeScale','on','HessUpdate','lbfgs',...
    'InitialHessType','identity','GoalsExactAchieve',1,...
    'TolFun',1e-10,'TolX',1e-9,'MaxIter',1000,...
'rho',.1,'sigma',0.900,'tau1',3,'tau2', 0.1,...
'tau3', 0.5,'StoreN',10);


objFunc = @(x) ObjGrad(x,X1,nrow, ncol,C,1);


x0 = x_opt;
%x0 = x(1:(nrow*ncol + n));
tic 
[x] = fminlbfgs(objFunc,x0,options);
toc

x = x';
lenLMI = nrow*ncol;
lenL = nrow/2;
% Get back required matrices:
U_opt = reshape(x(1:lenLMI,1),nrow, ncol) ;
L_opt = x(lenLMI+1:end,1);
V_opt = U_opt * U_opt';
P_opt = V_opt(1:n,1:n);
Q_opt = V_opt(1:n,n+1:end) - L_opt*C*P_opt;
A_opt = (Q_opt * inv(P_opt));
B_opt = X(:,n+1:end); 


% The output equation:
Rx = Rinv([X1_til;U1])*blkdiag(P_opt, 1); 
Y_i = Y*Rx*Rinv([X1_til;U1]*Rx);
C_opt = Y_i(1:n);
D_opt = Y_i(n+1:end);

x_id = zeros(n,length(t));
y_id = zeros(1,length(t)-1);
x_id(:,1) = 10*ones(n,1);
for i = 1:(length(t)-1)
    x_id(:,i+1) = A_opt*x_id(:,i) + B_opt;
    y_id(:,i) = C_opt*x_id(:,i) + D_opt;
end
figure(1);
errOut= 100*vecnorm((X_ts.' - x_id))./(vecnorm(X_ts.'));
plot(errOut)


figure(2);
plot(X_ts(:,600));
hold on
plot(x_id(600,:),'*');

