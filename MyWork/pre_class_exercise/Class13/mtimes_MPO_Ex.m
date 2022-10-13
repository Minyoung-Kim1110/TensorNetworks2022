function C = mtimes_MPO_Ex (B,A,Nkeep,Nsweep)
% < Description >
%
% C = mtimes_MPO_Ex (B,A,Nkeep,Nsweep)
%
% Variational multiplication of two matrix product operators (MPOs) B and
% A, where A multiplies to B from above (i.e., A is the top row, and B is
% the row below A). In the variational scheme, the multiplication result C
% is initialized as A, and then undergoes a two-site update, as described
% in App. D 1 of [B.-B. Chen et al., Phys. Rev. X 8, 031082 (2018)].
%
% < Input >
% B, A : [1 x L cell array] MPOs to be multiplied together. Each tensor
%       A{n} or B{n} is a rank-4 tensor acting on site n. Its legs are
%       ordered as bottom-top-left-right, where the bottom (top) leg
%       contracts to the physical leg of bra (ket) tensor. 
% Nkeep : [numeric] The maximum bond dimension for the result MPO C.
% Nsweep : [numeric] Number of round trips of sweeping. There will be
%       Nsweep pairs of left-to-right sweep and right-to-left sweep. That
%       is, the first sweep is left-to-right and the last is right-to-left.
%
% < Output >
% C : [1 x L cell array] Multiplication of A and B. Each tensor C{n}
%       follows the same leg order convention as thoes in A and B.
%
% Written by S.Lee (Oct.10,2022)

N = numel(A);

% sanity check
if N ~= numel(B)
    error('ERR: Length of two input MPOs do not match.');
end
for itN = (1:N)
    if (itN == 1) && ~all([size(A{itN},3),size(B{itN},3)] == 1)
        error('ERR: The leftmost leg of an MPO should be dummy.');
    elseif (itN == N) && ~all([size(A{itN},4),size(B{itN},4)] == 1)
        error('ERR: The rightmost leg of an MPO should be dummy.');
    elseif (itN < N) && (size(A{itN},4) ~= size(A{itN+1},3))
        error(['ERR: The fourth (= right) leg of A{',sprintf('%i',itN), ...
            '} and the third (= left) leg of A{',sprintf('%i',itN+1), ...
            '} do not have the same dimensions.']);
    elseif (itN < N) && (size(B{itN},4) ~= size(B{itN+1},3))
        error(['ERR: The fourth (= right) leg of B{',sprintf('%i',itN), ...
            '} and the third (= left) leg of B{',sprintf('%i',itN+1), ...
            '} do not have the same dimensions.']);
    elseif size(A{itN},1) ~= size(B{itN},2)
        error(['ERR: The first (= bottom) leg of A{',sprintf('%i',itN), ...
            '} and the second (= top) leg of B{',sprintf('%i',itN), ...
            '} do not have the same dimensions.']);
    end
end
% % % 

% Initialize C with A
C = A;

% Bring C into right-canonical form
% First, convert rank-4 tensors into rank-3, by merging physical legs, to
% use the canonForm function that canonicalize MPSs
Aloc = cell(1,N); % isometries for merging the bottom and top legs of MPO tensors
for itN = (1:N)
    Aloc{itN} = getIdentity(C{itN},1,C{itN},2);
    C{itN} = contract(C{itN},4,[1 2],Aloc{itN},3,[1 2]);
    
end
% Use canonForm for MPS
C = canonForm(C,0,Nkeep,[]);
% Bring back to rank-4
for itN = (1:N)
    C{itN} = contract(C{itN},3,3,conj(Aloc{itN}),3,3,[3 4 1 2]);
end

% Contractions of A, B, and C^\dagger. They correspond to the "effective
% Hamitonian" for the left and right parts in the ground and excited states
% search within the DMRG. As in the DMRG, we set ABC{1} and ABC{N+2} as
% 1's, to simplify the code. 
ABC = cell(1,N+2);
ABC{1} = 1;
ABC{end} = 1;
% % % % TODO (start) % % % %
% Feel free to define sub-functions for tensor contractions that appear
% multiple times in the code. By doing so, you will be able to write a
% simpler code with less probability of encountering bugs!

% calculate AB first 
%AB = cell(1, N);
%for itN = (1:N)
%    AB{itN} = contract(A{itN}, 4, 1, B{itN}, 4, 2, [4 1 2 5 3 6]);
%end 


% contract from right, since the first sweep is left-to-right
for itN = (N:-1:1)
    %ABC{itN+1}  = updateVR(A{itN}, B{itN}, C{itN}, ABC{itN+2});
    ABC{itN+1}  = updateVR(A{itN}, B{itN}, C{itN-1}, ABC{itN+2});
    
    %ABC{itN+1}  = updateVR(AB{itN},C{itN}, ABC{itN+2});
end

for itS = (1:Nsweep)
    isleft = true;
    for itN = (1:(N-1)) % left-to-right sweep      % change N-2 to N-1  
        [C{itN}, C{itN+1}] = update_2site(ABC{itN}, A{itN}, B{itN}, A{itN+1}, B{itN+1}, ABC{itN+3}, Nkeep, isleft);
        ABC{itN+1} = updateVL(A{itN}, B{itN}, C{itN}, ABC{itN});
        %[C{itN}, C{itN+1}] = update_2site(ABC{itN}, AB{itN}, AB{itN+1},  ABC{itN+3}, Nkeep, isleft);
        %ABC{itN+1} = updateVL(AB{itN}, C{itN}, ABC{itN});
        [C{itN}, C{itN+1}] = update_2site(ABC{itN}, ABC{itN+3}, Nkeep, isleft);
        ABC{itN+1} = updateVL(A{itN}, B{itN}, C{itN}, ABC{itN});
        
    end

    %   
    isleft = false;
    for itN = ((N-1):-1:1) % right-to-left sweep
        [C{itN}, C{itN+1}] = update_2site(ABC{itN}, A{itN}, B{itN}, A{itN+1}, B{itN+1}, ABC{itN+3}, Nkeep, isleft);
        ABC{itN+2}  = updateVR(A{itN+1}, B{itN+1}, C{itN+1}, ABC{itN+3});
        %[C{itN}, C{itN+1}] = update_2site(ABC{itN}, AB{itN},AB{itN+1},    ABC{itN+3}, Nkeep, isleft);
        %ABC{itN+2}  = updateVR(AB{itN+1}, C{itN+1}, ABC{itN+3});
    end
end
% % % % TODO (end) % % % %

end

function VR_new = updateVR(A, B, C, VR_old)
    VR_new  = contract(A, 4, 4, VR_old, 3, 1);
    VR_new = contract( VR_new, 5, [1 4],B, 4, [2 4]);
    C_dagger = permute(conj(C), [2 1 3 4]);
    VR_new = contract(VR_new, 5, [1 3 4], C_dagger, 4, [1 4 2]);
end 

function VL_new = updateVL(A, B, C, VL_old)
    VL_new  = contract(VL_old, 3, 1, A, 4, 3);
    VL_new = contract( VL_new, 5, [1 3], B, 4, [3 2]);
    C_dagger = permute(conj(C), [2 1 3 4]);
    VL_new = contract(VL_new, 5, [1 2 4], C_dagger, 4, [3 1 2]);
end 


function [left, right] = update_2site(VL, A1, B1, A2, B2, VR, Nkeep, isleft)
    new_state = contract(VL, 3, 1, A1, 4, 3);
    new_state= contract( new_state, 5, [1 3], B1, 4, [3 2], [1 4 2 3 5]);
    new_state = contract(new_state, 5, 4, A2, 4, 3);
    new_state = contract(new_state, 7, [4 5], B2, 4, [3 2]);
    new_state = contract(new_state, 7, [5 7], VR, 3, [1 2], [1 2 3 5 4 6]);
    [u, s, vd] = svdTr(new_state, 6, [1 2 3], Nkeep, []);
    if isleft
        left = permute(u, [2 3 1 4]);
        right = contract(diag(s), 2, 2, vd, 4, 1, [2 3 1 4]);
    else 
        left = contract(u, 4, 4, diag(s), 2, 1, [2 3 1 4]);
        right = permute(vd, [2 3 1 4]);
    end 
end 

%function VR_new = updateVR(AB, C, VR_old)
%    VR_new = contract(AB, 6, [5 6], VR_old, 3, [1 2]);
%    C_dagger = permute(conj(C), [2 1 3 4]);
%    VR_new = contract(VR_new, 5, [1 2 5], C_dagger, 4, [2 1 4]);
% end

%function [left, right] = update_2site(VL, cen1, cen2, VR, Nkeep, isleft)
%    new_state = contract(VL, 3, [1 2], cen1, 6, [3 4]);
%    new_state = contract(new_state, 5, [4 5], cen2, 6, [3 4]);
%    new_state = contract(new_state, 7, [6 7], VR, 3, [1 2]);
%    [u, s, vd] = svdTr(new_state, 6, [1 2 3], Nkeep, []);
%    if isleft
%        left = permute(u, [2 3 1 4]);
%        right = contract(diag(s), 2, 2, vd, 4, 1, [2 3 1 4]);
%    else 
%        left = contract(u, 4, 4, diag(s), 2, 1, [2 3 1 4]);
%        right = permute(vd, [2 3 1 4]);
%    end 
%end 

%function VL_new = updateVL(AB, C, VL_old)
%    VL_new  = contract(VL_old, 3, [1 2], AB, 6, [3 4]);
%    C_dagger = permute(conj(C), [2 1 3 4]);
%    VL_new = contract(VL_new, 5, [1 3 2], C_dagger, 4, [3 1 2]);
%end 
