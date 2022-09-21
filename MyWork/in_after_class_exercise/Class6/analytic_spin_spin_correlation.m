function [expectation_value] = analytic_spin_spin_correlation(alpha, beta, L, n)
%SPIN_SPIN_CORRELATION 이 함수의 요약 설명 위치
%   자세한 설명 위치

expectation_value= ( (-4/9)-4*( -1)^( alpha + beta) * (-1/3)^L);
expectation_value = expectation_value./(1+(-1)^(alpha+beta)*(-1/3)^L);

end

