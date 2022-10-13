L = 4; 
[S, I] = getLocalSpace('Spin', 1);

% bulk tensor for S_l \dot S_(l+1) 
%    | I                                                |
%    | Sp/sqrt(2)                                       |
%    | Sz                                               |
%    | Sm/sqrt(2)                                       |
%    | 0            Sp'/sqrt(2)    Sz'   Sm'/sqrt(2)  I |


Hloc = cell(14 ,14);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
for ito = (1:size(S,3)) % different components of spin operators
    Hloc{ito+1,1} = S(:,:,ito);
    Hloc{end,ito+1} = S(:,:,ito)';
end
for ito = (1:size(S,3))
    for ito2 = (1:size(S, 3))
        Hloc{ito + ito2-1,1} =S(:,:,ito) * S(:, :, ito2);

        Hloc{end, ito + ito2-1} =1/3*S(:,:,ito) * S(:, :, ito2);
    end 
end 

Hloc{end,end} = I;
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)])); % size = [2 2 5 5]

% bulk tensor for all 
Hs = cell(1,L);
Hs(:) = {Hloc};

Hs{1} = Hs{1}(:,:,end,:); % choose the last index of the left leg
Hs{end} = Hs{end}(:,:,:,1); % choose the first index of the right leg
