function alpha = backtrackingLineSearch(objFunc,objFuncValue,x,dx,dir)

%backtracking: this is also expensive for CFD but can be improved
%
% objFunc      - handle for objective function
% objFuncValue - current objective function value @ x
% x            - x
% dx           - dx
% dir          - search direction
%
% example : mb_backtrackingLineSearch(objFunc,objFuncValue,x,dx,dir)

alphaMax     = 1; % this is the maximum step length
alpha        = alphaMax;
fac          = 1/2; % < 1 reduction factor of alpha
c_1          = 1e-1;

while objFunc(x+alpha*dir) > objFuncValue + c_1*alpha*dir'*dx;

    alpha = fac*alpha;
 
end

end