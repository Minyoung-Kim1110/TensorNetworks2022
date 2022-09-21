function [expectation_value] = analytic_magnetization(alpha,beta, L, n)
%alpha, beta: boundary conditions 1,2 
% L : number of sites 
% N : site of magnetization 

expectation_value = 2 * (-1)^alpha * ((-1/3)^n - (-1)^(alpha + beta)*(-1/3)^(L-n+1))/(1+(-1)^(alpha + beta) * ( -1/3)^L);
end

