clear
L=6;
dt = 0.01;
[S, I] = getLocalSpace('Spin', 1/2);

% interaction term 
In =cat(3, S(:, :, 1),S(:, :, 3));

%Hamiltonian 
H = -contract(In, 3, 3, permute(conj(In), [2, 1, 3]), 3, 3);

% diagonalize and exponentiate 
ldim = size(H, 1);
Hmat = reshape(permute(H,[1 3 2 4]),(ldim^2)*[1 1]); % matrix representation
[V, D] = eig((Hmat + Hmat')/2);
D = diag(D);
expH = V * diag(exp(-1j * dt * D)) * V';
expH = reshape(expH, ldim*ones(1,4));

% two-site gate into two single-site tensor 
[U, d, V] = svdTr(expH, 4, [1, 3], [], []);

%left = contract(U, 3, 3, diag(d), 2, 1); % leg 1 and 2
%right = permute(V, [2 3 1]); % leg 3 and 4 
left = U; 
right = contract(diag(d), 2, 2, V, 3, 1, [2 3 1]);
% get MPO
Hs = cell(1,L);
Hs{1} = reshape(left, 2, 2, 1, []); 
for itN=(2:2:L-1)
    Hs{itN} = contract(right, 3, 2, left, 3, 1, [1 3 2 4]);
end
for itN=(3:2:L-1)
    Hs{itN}=contract(left, 3, 2, right, 3, 1, [1 3 4 2]);
end 
Hs{end} = reshape(right, 2, 2, size(right, 3), 1);

% check 
Hs_tot = 1; % initialize
for itN = (1:L)
    Hs_tot = contract(Hs_tot, 2*itN, 2*itN, Hs{itN},4,3);
end
Hs_tot = permute(Hs_tot,[(2:2*L+2) 1]);
Hs_tot = permute(Hs_tot,[(1:2:2*L) (2:2:2*L)]);
Hs_tot = reshape(Hs_tot,(size(I,1)^L)*[1 1]);

% construct exact Hamiltonian with kron product 
% just for L=6
expH_mat = reshape(expH, ldim*ldim, []);
odd = kron(kron(expH_mat, expH_mat), expH_mat);
even =kron(kron(kron(I, expH_mat), expH_mat), I);
Hamiltonian = odd * even;

% by taylor expansion 
%bulk tensor 
Hloc = cell(4,4);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = S(:,:,1);
Hloc{end,2} = 1j *dt* S(:, :, 1)';
Hloc{3,1} = S(:, :, 3);
Hloc{end,3} = 1j *dt* S(:, :, 3)';
Hloc{end,end} = I;
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)]));

% MPO for the full chain
Hs_t = cell(1,L);
Hs_t(:) ={Hloc};
Hs_t{1} = Hs_t{1}(:,:,end,:); % choose the last index of the left leg
Hs_t{1}(:, :, :, 1) = I;
Hs_t{end} = Hs_t{end}(:,:,:,1); % choose the first index of the right leg

% check 
Hs_tot_t = 1; % initialize
for itN = (1:L)
    Hs_tot_t = contract(Hs_tot_t, 2*itN, 2*itN, Hs_t{itN},4,3);
end
Hs_tot_t = permute(Hs_tot_t,[(2:2*L+2) 1]);
Hs_tot_t = permute(Hs_tot_t,[(1:2:2*L) (2:2:2*L)]);
Hs_tot_t = reshape(Hs_tot_t,(size(I,1)^L)*[1 1]);

max(max(Hamiltonian - Hs_tot)) 
max(max(Hamiltonian - Hs_tot_t))