J=1; 
L = 6; 
M = cell(1, L );
[S, I] = getLocalSpace('Spin', 1/2);

Hprev = 0 ; 
Aprev = 1; 

% For first site 
Anow = getIdentity(Aprev, 2, I, 2, [1,3,2]);
Hnow = updateLeft(Hprev, 2, Anow, [], [], Anow);
Sprev = updateLeft([],[],Anow, S, 3, Anow); % S operator in new basis 
Aprev = Anow; 
Hprev = Hnow;
% other sites 
for itN=(2:L)
    Anow = getIdentity(Aprev, 2, I, 2, [1,3,2]);
    Hnow = updateLeft(Hprev, 2, Anow, [], [], Anow);
%     spin-spin interaction  S^+S 
    Hsp = updateLeft(Sprev, 3, Anow, permute(conj(S), [2,1,3]), 3, Anow);
    Hnow = Hnow + J* Hsp ; 
    Aprev = Anow; 
    Hprev = Hnow; 
    if itN==L 
%         Diagonalize at last 
        [V, D] = eig((Hnow+Hnow')/2);
        [E_G, minid]  = min(diag(D));
        M{itN}=contract(Anow, 3, 2, V(:, minid), 2, 1, [1 3 2]);
    else 
        M{itN}=Anow; 
    end 

    Sprev = updateLeft([],[],Anow, S, 3, Anow); % S operator in new basis 
    disptime(['#', sprintf('%02i/%02i', [itN, L]), ' : NK=',sprintf('%i/%i', [size(M{itN}, 2), size(Hnow, 2)])]);
end 


% bulk tensor for each chain site
%    | I                                                |
%    | Sp/sqrt(2)                                       |
%    | Sz                                               |
%    | Sm/sqrt(2)                                       |
%    | 0           JSp'/sqrt(2)   JSz'  JSm'/sqrt(2)  I |


Hloc = cell(5,5);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
for ito = (1:size(S,3)) % different components of spin operators
    Hloc{ito+1,1} = S(:,:,ito);
    Hloc{end,ito+1} = J*S(:,:,ito)';
end
Hloc{end,end} = I;
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)])); % size = [2 2 5 5]

% bulk tensor for all 
Hs = cell(1,L);
Hs(:) = {Hloc};

Hs{1} = Hs{1}(:,:,end,:); % choose the last index of the left leg
Hs{end} = Hs{end}(:,:,:,1); % choose the first index of the right leg

% Get high rank taensor from MPO
Hs_tot = 1; % initialize
for itN = (1:L)
 Hs_tot = contract(Hs_tot, 2*itN , 2*itN ,Hs{itN}, 4,3);
 %   1          2*itN   3    |(2)   4
 % ---- Hs_tot ------ * ------Hs{itN}---
 %                           |(1)
end

% permute the left- and rightmost legs to the end
Hs_tot = permute(Hs_tot,[(2:2*L+2) 1]); % Same as (It implicitly erase last index since it is 1) Hs_tot = squeeze(Hs_tot);
% merge the incoming legs into a thick incoming leg;
% merge the outgoing legs into a thick outgoing leg
Hs_tot = permute(Hs_tot,[(1:2:2*L) (2:2:2*L)]);
Hs_tot = reshape(Hs_tot,(size(I,1)^L)*[1 1]); % 64 * 64 

HM = cell(1,L);
for itN = (1:L)
    %   1           2   
    % ---- M{itN} ------ 
    %       |
    %       3                     1              2
    %       |                   ======HM{itN}========
    %       *              =>            |
    %       |                            3
    %       2                            |
    %       |
    %   3   |      4   
    % ---- Hs{itN} ------ 
    %       |
    %       1
    %       |
    
    HM{itN} = contract(Hs{itN},4,2,M{itN},3,3);
    % leg order: Hbottom-Hleft-Hright-Mleft-Mright

    % isometry to merge left legs
    if itN == 1
        Aleft = 1; % there are only dummy legs
    else
        % use Aright from the previous iteration, 
        % to be a valid insertion of identity
        Aleft = conj(Aright);
    end
    
    % isometry to merge right legs
    % generate Identity
    %   1           2        2
    % ---- M{itN} ------ * --------
    %       |                      |
    %       3                      |
    %       |                      |    3
    %       *                    Aright------     
    %       |                      |      
    %       2                      |      
    %       |                      |
    %   3   |      4          1    |
    % ---- Hs{itN} ------ * -------
    %       |
    %       1
    %       |
    Aright = getIdentity(HM{itN},3,HM{itN},5);
    
    % contract isometries
    %         2     1           2     2 
    %        ----*---- M{itN} ----* ----  
    %       |            |               |
    %       |            3               |                  
    %       |            |               |        
    %   3   |            *               |      3         
    %---- Aleft          |             Aright ----                  
    %       |            2               |                 
    %       |            |               |
    %       |  1     3   |      4     1  | 
    %        ----*---- Hs{itN} ----* ---- 
    %                    |
    %                    1
    %                    |
    HM{itN} = contract(Aleft,3,[1 2],HM{itN},5,[2 4]);
    % leg order: Aleft-Hbottom-Hright-Mright
    HM{itN} = contract(HM{itN},4,[3 4],Aright,3,[1 2],[1 3 2]);
end
[HM,HMnorm] = canonForm(HM,L,[],[]) % left-canonical
[HM,HMnorm2] = canonForm(HM,0,[],[]) % here the second output should be 1
HM{1} = HM{1}*HMnorm;
MHM = 1;
for itN = (1:L)
    MHM = updateLeft(MHM,2,M{itN},[],[],HM{itN});
end
MHM - E_G % zero up to numerical noise