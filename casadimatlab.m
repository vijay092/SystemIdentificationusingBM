
addpath('C:\Users\sanja\Google Drive\casadi-windows-matlabR2016a-v3.5.1')
import casadi.*
clear all;
close all;
set(0,'DefaultFigureWindowStyle','docked')


n = 100;
ts = 0.5;
csys = rss(n);
sys = c2d(csys, ts);


% Full order system matrices: 
t = 0:ts:10;
X = zeros(n,length(t));
Y = zeros(1,length(t)-1);
A = sys.A; B = sys.B; C = sys.C; D = sys.D;
for i = 1:(length(t)-1)
    X(:,i+1) = A*X(:,i) + B;
    Y(:,i) = C*X(:,i) + D;  
end 
figure(1);
plot(Y);

% Construct snapshots:
X1_til = X(:,1:end-1);
X2_til = X(:,2:end);
u = ones(length(t),1)';
U1 = u(1:end-1);

% Compute R inverse:
X = X2_til*Rinv([X1_til;U1]);
X1 = X(1:n,1:n);

% Weighting matrices:
X1 = X1;
X2 = [eye(n); zeros(n)];
X3 =  X2';
X4 = [zeros(n);eye(n)];

% Optimization
nrow = 2*n;
ncol = 2*n;
x0 = 0.1*rand(nrow*ncol,1);
fun = @(x) CostFn(x,X1,X2, X3, X4, nrow,ncol);
options = optimoptions(@fminunc,'MaxFunctionEvaluations',4e5,...
            'MaxIterations',1e3, 'PlotFcn', 'optimplotfval',...
            'OptimalityTolerance',1e-10,...
            'StepTolerance',1e-10)
tic        
x_opt = fminunc(fun,x0,options);
toc
% Get back required matrices:
U_opt = reshape(x_opt,nrow,ncol);
V_opt = U_opt * U_opt';
P_opt = V_opt(1:n,1:n);
Q_opt = V_opt(1:n,n+1:end);
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
figure(1)
hold on;
plot(y_id,'r*')




