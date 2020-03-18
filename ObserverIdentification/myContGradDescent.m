% =========================================================================
%% This code identifies an ob




clear all; close all
n = 5;
csys = rss(n); 
ts = .5;
sys = c2d(csys, ts);


% Full order system matrices: 
t = 0:ts:10;
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

% Compute R inverse:
X = X2_til*Rinv([X1_til;U1]);
X1 = X(1:n,1:n);
nrow = 2*n;
ncol = 2*n;

% In this section, I have computed the analytical gradient and instead of
% doing the vanilla gradient descent, I used the continuous version of it
% because i was having problems with constant step size-- I am working on
% adaptive step size now. 
gradFun = @(t,x) gradODE(t,x,X1,nrow,ncol,C);


%options = odeset('RelTol',1e-6,'AbsTol',1e-8);
x0 = rand(nrow*ncol + nrow/2,1)
[t,x_opt] = ode23(gradFun,t,x0);

% Get back required matrices:
U_opt = reshape(x_opt(end,1:nrow*ncol),nrow, ncol) ;
L_opt = x_opt(end, nrow*ncol+1:end)
V_opt = U_opt * U_opt';
P_opt = V_opt(1:n,1:n);
Q_opt = V_opt(1:n,n+1:end) - L_opt'*C*P_opt;
A_opt = (Q_opt * inv(P_opt));
B_opt = X(:,n+1:end); 


% The output equation:6
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



