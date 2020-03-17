% This code obtains an analytical gradient for the cost that I am
% interested in


addpath('C:\Users\sanja\Google Drive\casadi-windows-matlabR2016a-v3.5.1')
import casadi.*
clear all;
close all;
set(0,'DefaultFigureWindowStyle','docked')


% n = 1;
% ts = 0.5;
% csys = rss(n);
% sys = c2d(csys, ts);
% 
% 
% % Full order system matrices: 
% t = 0:ts:20;
% X = zeros(n,length(t));
% Y = zeros(1,length(t)-1);
% A = sys.A; B = sys.B; C = sys.C; D = sys.D;
% for i = 1:(length(t)-1)
%     X(:,i+1) = A*X(:,i) + B;
%     Y(:,i) = C*X(:,i) + D;  
% end 
% figure(1);
% plot(Y);
% 
% % Construct snapshots:
% X1_til = X(:,1:end-1);
% X2_til = X(:,2:end);
% u = ones(length(t),1)';
% U1 = u(1:end-1);
% 
% % Compute R inverse:
% X = X2_til*Rinv([X1_til;U1]);
% X1 = X(1:n,1:n);
% 
% % Weighting matrices:
% X1 = [X1 zeros(n)];
% X2 = [eye(n); zeros(n)];
% X3 =  X2';
% X4 = [zeros(n);eye(n)];

% syms U [2*n,2*n] real
% 
% grad = gradient((norm(X1*U*X2 - X3*U*X4 ,'fro'))^2,reshape(U,2*n*2*n,1));
% grad_zero =  solve(grad == 0);
% grad_zero.U1_1


syms P1 [2,2] real
syms Q1 [2,2] real
syms P2 [2,2] real
syms Q2 [2,2] real
syms dell [2,2] real

V = reshape([P1, Q1; Q2 P2],16,1)
reshape(jacobian(gradient(norm(P1-Q1,'fro')^2,V),V),16,16)





