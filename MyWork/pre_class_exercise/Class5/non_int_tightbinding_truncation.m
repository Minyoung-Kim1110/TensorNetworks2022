clear
% number of sites 
% number of energy states to keep 
Nkeep=300;
Nmax = 50 ; 
result = zeros(4, Nmax);

% (i)
for N=(2:Nmax)
    t = ones(1, N-1);
    result(1, N) =non_intTB_truncation(t, Nkeep);
    result(2, N) = nonIntTB(t);
end

% (ii)
for N=(2:Nmax)
    t = 2.^((0:N-1)/2);
    result(3, N) = non_intTB_truncation(t, Nkeep);
    result(4, N) = nonIntTB(t);
end
result