
function A = Wstate(N)
if N == 1
    A = 1; 
elseif N==2 
    A = 1/sqrt(2) * (ones(2) - eye(2));
else
    A = ones(2) - eye(2);
    for i=3:N
        s = size(A);
        dim = prod(s, 'all');
        add = zeros(1, dim);
        add(1) = 1; 
        add = reshape(add, s);
        A = cat(i, A, add);
    end 
    A = A /sqrt(N); 
end 