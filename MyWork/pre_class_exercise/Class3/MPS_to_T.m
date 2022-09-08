function tensor = MPS_to_T(Q)
%From matrix product states, contract to a high-rank tensor. 
% < input >
%    Q : [cell(n)] each component contains matrix product states 
%    
% < output >
%    tensor: [numeric tensor] a high-rank tensor contracted from Q(rank n) 
% Written by M.Kim(Sep.07,2022)

n = numel(Q);
tensor = Q{1};
for it = (2:n)
    rank = numel(size(tensor));
    tensor = contract(tensor, rank, rank-1, Q{it}, 3, 1);
end
shape = size(tensor);
shape(numel(shape)-1) = shape(numel(shape));
shape(numel(shape)) = shape(numel(shape)-1);
tensor = reshape(tensor, shape(2:end-1)); % Drop first dimension, last dimension since it is dummy
end

