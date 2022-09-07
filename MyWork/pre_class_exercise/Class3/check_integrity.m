function same = check_integrity(tensor)
%From tensor, decompose to matrix product states and contract the states.
%Check whether original tensor is reconstructed.
% < input >
%    tensor: [numeric tensor] high-rank tensor to be decomposed(rank n) 
% < output >
%    same: [integer/boolean] 1 if original tensor is constructed else 0
% WrittTen by M.Kim(Sep.07,2022)


if abs(sum(tensor, 'all')-1)>10^-8
    tensor = tensor/norm(tensor(:));

Q = T_to_MPS(tensor);
T = MPS_to_T(Q);

num_element = numel(T);
same = abs(T-tensor)<10^-8;
same = (num_element == sum(same, 'all'));

end

