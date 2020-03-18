% Using the closed form solution of gradient:
% addpath('C:\Users\sanja\Google Drive\casadi-windows-matlabR2016a-v3.5.1')
% clear all;
% close all;
% set(0,'DefaultFigureWindowStyle','docked')
% 
% 
% n = 10;

% csys = rss(n);
% sys = c2d(csys, ts);
clear all; close all
load stateMat

% Full order system matrices: 
n = length(A)
D = -0.1;
ts = .5;
t = 0:ts:40;
X = zeros(n,length(t));
Y = zeros(1,length(t)-1);
%A = sys.A; B = sys.B; C = sys.C; D = sys.D;
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

% Gradient Descent:
alpha = .1;
U = eye(2*n);
j = 0; eps = 10;
while eps > 5e-1
                                                                                                        
    V = U*U';
    P = V(1:n,1:n);
    Q = V(1:n,n+1:end);
    grad = 4*[X1'*X1*P - Q'*X1, zeros(n) ;...
               (-X1*P + Q)  , zeros(n)]*U; 
    U = U - alpha * grad;
    eps = norm(grad,'fro'); 
    j = j +1;
    fprintf('It = %f , \t Grad = %f\n',j,eps) 
                 
    if j > 700000
        break;
    end 
    alpha = 0.7*alpha;                                                            

end 

% Get back required matrices:
U_opt = U;
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


