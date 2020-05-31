
clear all
load StateMatricesDoubleVortex
ns = 2380;

A = A(1:ns,1:ns);
b = b(1:ns);
n = length(A);
ny = 1;
% Generate time-series data
t = 0:1:500;
x0 = 0*ones(n,1);
[t,X_ts] = ode45(@(t,x) A*x + b, t, x0);

plot(t, X_ts);
xtrim = X_ts(end,:);



idx = find(xtrim> 200);



save('xtrim','xtrim','idx')