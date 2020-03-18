function R = Rinv(V)

R = V'*pinv(V*V');

end