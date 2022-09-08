function [M, entropys] = T_SVD(tensor)

sz = size(tensor);
M = cell(1,numel(sz));
entropys = cell(1, numel(sz)-1);
A = tensor; 
remain = prod(sz);
szl = 1; % the bond dimension of the left leg of Q{n}for it = (1:(numel(sz)-1))
for it = (1:numel(sz)-1)
    remain = remain/sz(it);
    A = reshape(A,[szl*sz(it), remain]);
    [U,S,V] = svd(A, 0);
    s = diag(S);
    s=s(abs(s)>10^(-16));
    Spart = (s.*s).* log2(s.*s);
    entropys{it} = - sum(Spart(~isnan(Spart)));
    M{it} = permute(reshape(U, szl, sz(it), []), [1,3,2]);
    szl = size(M{it}, 2);
    A = reshape(S* V' , [szl, sz(it+1:end)]);
    

end
M{end} = permute(A, [1,3,2]);
end
