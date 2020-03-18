clc, clear; close all;
rng(1231);
OPTIM_CHC = 1; % 0: gradient descent without momentum, 1: with momentum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: S. Vijayshankar, A. Chakrabarty
% Creation Date: Mar 6, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 50; % state space dim
csys = rss(n); % random state-space CT model
ts = 1;
sys = c2d(csys, ts); % convert to DT

%% Forward simulations to get time series data
% Snapshots: [X: states, Y: outputs]
t = 0:ts:20;
X = zeros(n,length(t));
Y = zeros(1,length(t)-1);
A = sys.A; B = sys.B; C = sys.C; D = sys.D;
n = length(A);
for i = 1:(length(t)-1)
    X(:,i+1) = A*X(:,i) + B;
    Y(:,i) = C*X(:,i) + D;  
end 

% Construct snapshots:
X1_til = X(:,1:end-1);
X2_til = X(:,2:end);
u = ones(length(t),1)';
U1 = u(1:end-1);

% Compute R inverse:
X = X2_til*Rinv([X1_til;U1]);
X1 = X(1:n,1:n);
X12 = X1.'*X1;
nX1 = numel(X1);

%% Optimization [Adam - adv. ver. of SGD for CNN training]

U = randn((2*n)^2,1); % initial guess


if OPTIM_CHC == 1
    v = 0; gamma = 0.9;
end
step_size = 0.01; % needs tuning based on number of dimensions
optim_iter_max = 1e6;
grad_tol = 4e-4;

grad = eval_gradient(U, X1, X12);
grad_norm = norm(grad);

for k = 1:optim_iter_max
    
    grad = (1/nX1)*eval_gradient(U, X1, X12);
    grad_norm = norm(grad);
    if grad_norm < grad_tol
        break;
    end
    
    if OPTIM_CHC == 0
        % vanilla gradient descent
        U = U - step_size*grad(:);
    else
        % gradient descent with momentum
        v = gamma*v + step_size*grad(:);
        U = U - v;
    end
    
    if mod(k, 500) == 1
        fprintf(1, 'Iter: %d, Gradient Norm: %.4e \n', k, grad_norm);
    end
    
    %{
    if k == 1e5
        step_size = step_size/2;
    elseif k == 1e6
        step_size = step_size/2;
    end
    %}
    
end

%% Validation
x_opt = U.';

% Get back required matrices:
U_opt = reshape(x_opt(end,:),2*n, 2*n);
V_opt = U_opt * U_opt';
P_opt = V_opt(1:n,1:n);
Q_opt = V_opt(1:n,n+1:end);
A_opt = Q_opt / P_opt;
B_opt = X(:,n+1:end); 

% The output equation:
Rx = Rinv([X1_til;U1])*blkdiag(P_opt, 1); 
Y2 = Y*Rx*Rinv([X1_til;U1]*Rx);
C_opt = Y2(1:n);
D_opt = Y2(n+1:end);
x_id = zeros(n,length(t));
y_id = zeros(1,length(t)-1);

for i = 1:(length(t)-1)
    x_id(:,i+1) = A_opt*x_id(:,i) + B_opt;
    y_id(:,i) = C_opt*x_id(:,i) + D_opt;
end
figure(1)
t = t(1:end-1);
plot(t, Y, 'k-o', t, y_id, 'r-', 'linewidth', 2);
legend('Measured', 'Estimated');

function grad = eval_gradient(U, X1, X12)

% shape
nu = sqrt(length(U)); n = nu/2;

U = reshape(U,nu,nu);
V = U*U';
P = V(1:n,1:n);
Q = V(1:n,n+1:end);

temp = -X1*P + Q;
grad = [(X12*P - Q.'*X1), temp ;...
            temp', (X12*P - Q'*X1)]*U; 

end 