clear
% Initialization 
J = 1; 
N = 3; 
[S, I] = getLocalSpace('Spin', 1/2);
Ss = cell(1, N);
H = 0;
Aprev =1;  

% for one particle 
itN=1;
Anow = getIdentity(Aprev, 2, I, 2, [1,3,2]);
H = updateLeft(H, 2, Anow, [], [], Anow);
%below for loop do nothing since there is only one particle itN =1! 
for itN2 =(1: itN-1)
    Hsp = updateLeft(Ss{itN2}, 3, Anow, permute(conj(S), [2  1 3]), 3, Anow);
    H = H + J * Hsp ;
end
for itN2 = (1:itN)
    if itN2<itN 
        Ss{itN2} = updateLeft(Ss{itN2}, 3, Anow, [], [], Anow);
    else
        Ss{itN2}=updateLeft([], [], Anow, S, 3, Anow);
    end
end
Aprev = Anow;

% for two particles 
itN=1;
% span Hilbert space 
Anow = getIdentity(Aprev, 2, I, 2, [1,3,2]);
% convert H with new basis of Hilbert space 
H = updateLeft(H, 2, Anow, [], [], Anow);

% HERE!!! 

% for itN2 =(1: itN-1)
    % Hsp = updateLeft(Ss{itN2}, 3, Anow, permute(conj(S), [2  1 3]), 3, Anow); 
    %H = H + J * Hsp ;
%end
% current spin oerator in new basis
SA = contract(Ss{1}, 3, 2, Anow, 3, 1, [1,3,2,4]);
size(SA) % it should be (2, 4, 3, 2)
%   2    2   4
%  -- S -- A--
%     | 3  |2  

A_SA = contract( permute(conj(Anow), [2, 1, 3]), 3, 2, SA, 4, 1, [1,3,2,4,5]); 
%   4   2    2   4
%  --A'-- S -- A--
%    | 2  | 3  |2  
size(A_SA) % it should be (4, 4, 2, 3, 2)
total =  contract(A_SA, 5, [3,4,5], permute(conj(S), [2,1,3]), 3, [1, 3, 2]); 
%   4   2    2   4
%  --A'-- S -- A--
%    | 2  | 3  |2  
%    |    |    |
%    -----S'---- 
size(total) % it should be (4, 4)



%for itN2 = (1:itN)
%    if itN2<itN 
%        Ss{itN2} = updateLeft(Ss{itN2}, 3, Anow, [], [], Anow);
%    else
%        Ss{itN2}=updateLeft([], [], Anow, S, 3, Anow);
%    end
%end
%Aprev = Anow;