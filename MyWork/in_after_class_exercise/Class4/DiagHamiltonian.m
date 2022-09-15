J = 1; 
sigmas = cell(1, 3);
% Sx
sigmas{1} = zeros(2);
sigmas{1}(1, 2) =1/2;
sigmas{1}(2, 1) =1/2;
% Sy
sigmas{2} = zeros(2);
sigmas{2}(1, 2) = -1i/2;
sigmas{2}(2, 1) = 1i/2;

% S_z 
sigmas{3} = zeros(2);
sigmas{3}(1,1) = 1/2; 
sigmas{3}(2,2) = -1/2;

H = zeros(2^3);
% generate Hamiltonian with kronecker Tensor product 
for it = (1:3)
    H = H+kron(kron(sigmas{it}, sigmas{it}), eye(2))+ kron(kron(sigmas{it}, eye(2)), sigmas{it})+ kron(kron(eye(2), sigmas{it}), sigmas{it});
end 
H = J * H; 

[V, diag] = eig(H);
diag




