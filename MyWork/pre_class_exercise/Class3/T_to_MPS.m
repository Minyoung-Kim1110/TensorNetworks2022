function Q = T_to_MPS(tensor)
% Decompose a high-rank tensor into a matrix product state that consists of
% rank -2 and -3 tensors using QR decomposition 
% < input >
%    tensor: [numeric tensor] high-rank tensor to be decomposed(rank n) 
% < output >
%    Q : [cell(n)] each component contains matrix product states 
% Written by M.Kim(Sep.07,2022)
% Adapted from TensorDecomposition.m by S.Lee

sz = size(tensor);
Q = cell(1,numel(sz));
R = tensor; 
remain = prod(sz);
szl = 1; % the bond dimension of the left leg of Q{n}for it = (1:(numel(sz)-1))
for it = (1:numel(sz)-1)
    remain = remain/sz(it);
    R = reshape(R,[szl*sz(it), remain]);
    [Q{it},R] = qr(R,0); %QR decomposition efficiently
    Q{it} = reshape(Q{it}, szl, sz(it), []);
    Q{it} = permute(Q{it}, [1,3,2]); % permute to the left-right-bottom order
    szl = size(Q{it},2); % update the bond dimension
    R = reshape(R, [szl,sz(it+1:end)]);
end
Q{end} = permute(R, [1,3,2]);
end

