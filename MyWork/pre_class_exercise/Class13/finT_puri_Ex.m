function [taus,lnZs,PTS] = finT_puri_Ex (Hs,dt,tmax,Nkeep,Nsweep)
% < Description >
%
% [taus,lnZs,PTS] = finT_puri_Ex (Hs,dt,tmax,Nkeep,Nsweep)
%
% This function obtains the purified representation of unnormalized thermal
% density matrix exp(-beta*H), which is not divided by the partition
% function Z = Tr[ exp(-beta*H)], at temperature 1/beta = 1/(2*tmax), by
% simulating the imaginary-time evolution of the purified state, starting
% from one at infinite temperatures (beta = 0).
% In the way of the imaginary-time evolution, this function also computes
% the logarithm of the partition functions at temperatures 1/2/taus,
% corresponding to individual time steps.
% Here we decompose the time evolution operator by using the second-order
% Trotter decomposition for time steps of size dt. At each iteration, it
% applies the matrix product operator (MPO) representation of the time
% evolution operator (for time step dt) to the MPO representation of the
% purified thermal state, by using the variational MPO multiplication
% implemented in "DMRG/mtimes_MPO.m".
%
% < Input >
% Hs : [cell] Hamiltonian. Each cell element Hs{n} describes the two-site
%       interaction between site n and n+1. Thus, Hs(1:2:end) acts on odd
%       bonds; Hs(2:2:end) on even bonds. numel(Hs) is the system size
%       minus 1.
%       The leg convention of Hs{n} are as follows:
%
%       2      4      [legs 1 and 2 are for site n;
%       |      |       legs 3 and 4 are for site n+1]
%      [ Hs{n}  ]
%       |      |
%       1      3
%
% dt : [numeric] Time step size. Each imaginary-time evolution step (for
%       the imaginary time step -1i*dt) is decomposed into three exponential
%       terms, exp(-dt/2*Hodd) * exp(-dt*Heven) * exp(-dt/2*Hodd), where
%       Hodd (Heven) is the sum of Hs(1:2:end) [Hs(2:2:end)].
% tmax : [numeric] Maximum time to be reached at the end of the time
%       evolution; the whole time evolution simulates exp(-tmax*H), where H
%       is the sum of Hs(:). Since we are using the purification method,
%       the value of "tmax" should be set as a half of the inverse
%       temperature.
% Nkeep : [integer] Maximum bond dimension of the MPO form of the purified
%       thermal state.
% Nsweep : [integer] Number of round trips in the variational MPO
%       multiplication. "Nkeep" and "Nsweep" are directly forwarded to
%       "DMRG/mtimes_MPO.m"; see the description of those inputs therein.
%
% < Output >
% taus : [numeric] Time instances taken by the Trotterized imaginary-time
%       evolution steps. The actual imaginary time instances are -1i*taus.
%       The corresponding values of the inverse temperature are given by
%       2*taus.
% lnZs : [numeric] Logarithms of the partition function measured at time
%       instances "taus".
% PTS : [cell] The MPO representation of the purified thermal state at the
%       last time instance.
%
% Written by S.Lee (Oct.10,2022)

tobj = tic2;

N = numel(Hs)+1;

% % % check the integrity of input
for itN = (1:numel(Hs))
    if size(Hs{itN},1) ~= size(Hs{itN},2) 
        error(['ERR: The first and second legs of Hs{', ...
            sprintf('%i',itN),'} should have the same dimensions.']);
    elseif size(Hs{itN},3) ~= size(Hs{itN},4) 
        error(['ERR: The third and fourth legs of Hs{', ...
            sprintf('%i',itN),'} should have the same dimensions.']);
    elseif (itN < numel(Hs)) && size(Hs{itN},3) ~= size(Hs{itN+1},1)
        error(['ERR: The third leg of Hs{', ...
            sprintf('%i',itN),'} and the first leg of Hs{',sprintf('%i',itN+1), ...
            '} should have the same dimensions.']);
    end
end
% % % 

% number of time steps
Nstep = ceil(tmax/dt);
taus = dt*(1:Nstep);
lnZs = zeros(size(taus)); % result

% show information
fprintf('Finite T: Purification\n');
fprintf(['N = ',sprintf('%i',N),', Nkeep = ',sprintf('%i',Nkeep), ...
    ', Nsweep = ',sprintf('%i',Nsweep), ...
    ', dt = ',sprintf('%.4g',dt),', tmax = ',sprintf('%g',taus(end)), ...
    ' (',sprintf('%.4g',Nstep),' steps)\n']);

% Obtain two-site time evolution operators from Hs{..}, and decompose them
% into pairs of one-site tensors.
Os = cell(3,numel(Hs)); % array of one-site tensors; 1st and 3rd rows are 
% from the two-site gates from Hodd; 2nd row from Heven.

for itN = (1:numel(Hs))
    % reshape nearest-neighbor interaction term and diagonalize it
    ds = [size(Hs{itN}),ones(1,4-ndims(Hs{itN}))];
    H_mat = reshape(permute(Hs{itN},[1 3 2 4]),[prod(ds([1 3])), prod(ds([2 4]))]);
    [V,D] = eig((H_mat+H_mat')/2);
    [D,ids] = sort(diag(D),'ascend');
    V = V(:,ids);

    % % % % TODO (start) % % % %
    % two-site gate for imaginary-time evolution; full step of dt for even,
    % half step of dt/2 for odd
    if mod(itN, 2) == 1 % odd 
        expH = V * diag(exp(- dt/2 *D))*V';
    else 
        expH = V * diag(exp(- dt * D))*V';
    end
    % decompose Udt into Ldt and Rdt
    expH = reshape(expH, ds);
    [u, s, vd] = svdTr(expH, 4, [1 3], [], []);
    Ldt = contract(u, 3, 3, diag(s), 2, 1);
    Rdt = permute(vd, [2 3 1]);
    % assign Ldt and Rdt to Os
    if mod(itN,2) == 1 % odd 
        Os{1, itN} = Ldt; 
        Os{1, itN+1 } = Rdt; 
        Os{3, itN} = Ldt; 
        Os{3, itN+1} = Rdt;
    else % even
        Os{2, itN} = Ldt; 
        Os{2, itN+1} = Rdt; 

    end
    % % % % TODO (end) % % % %
end

% fill in identities
Os{2,1} = getIdentity(Os{3,1},2);
if isempty(Os{1,end})
    if ~isempty(Os{3,end}) || isempty(Os{2,end})
        error('ERR: It should not happen!');
    end
    Os{1,end} = getIdentity(Os{2,end},2);
    Os{3,end} = getIdentity(Os{2,end},1);
end
if isempty(Os{2,end})
    if isempty(Os{1,end}) || isempty(Os{3,end})
        error('ERR: It should not happen!');
    end
    Os{2,end} = getIdentity(Os{3,end},2);
end

% % % % TODO (start) % % % %
% Contract tensors along the same column, to obtain the MPO representation
% for a time evolution step over dt
U_MPO = cell(1,N);
for itN = (1:N)
    if mod(itN,2) == 1
        % merge the right legs of the tensors Os([1 3],itN)
        U_MPO{itN} =  contract(Os{3, itN}, 3, 2, Os{2, itN}, 3, 1);
        U_MPO{itN} = contract(U_MPO{itN}, 4, 3, Os{1, itN}, 3, 1, [1 4 3 2 5]);
        U_MPO{itN} = reshape(U_MPO{itN}, size(U_MPO{itN}, 1),size(U_MPO{itN}, 2),size(U_MPO{itN}, 3), []);
    else
        % take the "Aright" generated from the last iteration, to merge the
        % left legs
        U_MPO{itN} =  contract(Os{1, itN}, 3, 2, Os{2, itN}, 3, 1);
        U_MPO{itN} = contract(U_MPO{itN}, 4, 3, Os{3, itN}, 3, 1, [1 4 2 5 3 ]);
        U_MPO{itN} = reshape(U_MPO{itN}, size(U_MPO{itN}, 1),size(U_MPO{itN}, 2),[], size(U_MPO{itN}, 5));
    end
end
% % % % TODO (end) % % % %

% initialize the purified state for thermal state at infinite temperature
PTS = cell(1,N);
for itN = (1:N-1)
    PTS{itN} = getIdentity(Hs{itN},2);
end
PTS{end} = getIdentity(Hs{end},4);
% NOTE: We consider unnormalized thermal state, to compute the partition
% function as the trace of the purified state times its Hermitian
% conjugate (= the contraction of the purified state and its Hermitian
% conjugate).

% normalize 
norm = get_Trace(PTS);
PTS{1} =PTS{1} /sqrt(norm); 

% perform imaginary-time evolution
for it1 = (1:Nstep)
    % apply the time evolution operator U_MPO to the purified state MPO
    PTS = mtimes_MPO_Ex(U_MPO,PTS,Nkeep,Nsweep);

    % % % % TODO (start) % % % %
    % trace of the purified state times its Hermitian conjugate (= the
    % contraction of the purified state and its Hermitian conjugate)
    
    z = get_Trace(PTS); % square norm of purified state 
    % add the log of the trace to the result array lnZs
    lnZs(it1) = log(z);  % Natural Logarithm
    % add the result from the last iteration, since the MPO gets normalized
    % at every iteration (see below) to avoid divergence    
    if it1>1
        lnZs(it1) = lnZs(it1) +lnZs(it1-1);
    else 
        lnZs(it1) =  lnZs(it1) +log(norm);
    end
    % normalize the MPO
    PTS{1} = PTS{1}/sqrt(z);

    % % % % TODO (end) % % % %

    if (mod(it1,round(Nstep/10)) == 0) || (it1 == Nstep)
        disptime(['#',sprintf('%i/%i',[it1,Nstep]), ...
            ' : t = ',sprintf('%g/%g',[taus(it1),taus(end)])]);
    end
end

toc2(tobj,'-v');
chkmem;

end

function ret = get_Trace(purified_state)
    N = numel(purified_state);
    itN = N;
    ret = contract(purified_state{itN}, 4, [1 2 4], conj(purified_state{itN}), 4, [1 2 4], [2 1]);
    for itN=(N-1:-1:1)
        tmp = contract(purified_state{itN}, 4, [1 2], conj(purified_state{itN}), 4, [1 2], [1 3 2 4]);
        ret= contract(tmp, 4, [3 4], ret, 2, [1 2]);
    end 
end 

%function  ret= get_density_matrix(purified_state)
%    N = numel(purified_state);
%    ret = cell(1, N);
%    for itN=(N:-1:1)
%        if itN ==N 
%            ret{itN} = contract(purified_state{itN}, 4, [1 4], conj(purified_state{itN}), 4, [1 4], [3 2 1 4]);
%            ret{itN} = reshape(ret{itN}, size(ret{itN}, 1),size(ret{itN}, 2), []);
%        else 
%            ret{itN} = contract(purified_state{itN}, 4, 1, conj(purified_state{itN}), 4, 1, [4 1 2 5 3 6]);
%            ret{itN} = reshape(ret{itN}, size(ret{itN}, 1), size(ret{itN}, 2),size(ret{itN}, 3)*size(ret{itN}, 4), []);
%        end 
%    end 
%end 


