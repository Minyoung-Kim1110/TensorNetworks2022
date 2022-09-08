% generate a rank 5 tensor T 
sz = [2 3 2 3 4];
T = reshape((1:prod(sz)), sz);
T = T/norm(T(:));

% reshape tensor to matrix A & B 

% system A = {1, 2} system B = {3,4,5}
T_a = reshape(T, sz(1)*sz(2), []);

% system A = {1, 3} system B = {2,4,5}
T_b = reshape(permute(T, [1,3,2,4,5]), sz(1)*sz(3), []);

% system A = {1, 5} system B = {2,3,4}
T_c = reshape(permute(T, [1,5,2,3,4]), sz(1)*sz(5), []);

% Use SVD and calculate entanglement entropy using singular values
% Note that singular values are non-negative real values 
s_a = svd(T_a);
entropy_a = - sum((s_a .*s_a) .* log2( (s_a .* s_a)));
s_b = svd(T_b);
entropy_b = - sum((s_b .*s_b) .* log2( (s_b .* s_b)));
s_c = svd(T_c);
entropy_c = - sum((s_c .*s_c) .* log2( (s_c .* s_c)));
