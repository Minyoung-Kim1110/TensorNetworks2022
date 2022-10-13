clear

% system parameter
J1 = 1; % coupling strength
J2 = 1/2; % next-nearest coupling strength

% DMRG parameter
Nkeep = 30; % bond dimension
Nsweep = 4; % number of pairs of left+right sweeps

% Local operators
[S,I] = getLocalSpace('Spin',1/2);

% % MPO formulation of Hamiltonian
% Hamiltonian tensor for each chain site
Hloc = cell(8,8);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = S(:,:,1);
Hloc{3,1} = S(:,:,2);
Hloc{4,1} = S(:,:,3);

Hloc{end-1, 4} = I;
Hloc{end-2, 3} = I;
Hloc{end-3, 2} = I;

Hloc{end,2} = J1*S(:,:,1)';
Hloc{end,3} = J1*S(:,:,2)';
Hloc{end,4} = J1*S(:,:,3)';
Hloc{end,5} = J2*S(:,:,1)';
Hloc{end,6} = J2*S(:,:,2)';
Hloc{end,7} = J2*S(:,:,3)';

Hloc{end,end} = I;
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)]));

L = 40; % number of sites in a chain


% full chain
Hs = cell(1,L);
Hs(:) = {Hloc};
Hs{1} = Hs{1}(:,:,end,:); % choose the last components of the left leg
Hs{end} = Hs{end}(:,:,:,1); % choose the first components of the right leg

Minit = cell(1,L);
Minit{1} = rand(1,Nkeep,size(I,2));
Minit{end} = rand(Nkeep,1,size(I,2));
for itN = (2:L-1)
    Minit{itN} = rand(Nkeep,Nkeep,size(I,2));
end

[M0_2site,E0_2site,Eiter_2site] = DMRG_GS_2site_Ex(Minit,Hs,Nkeep,Nsweep);
E0_exact = - 3*L/8;
plot((1:numel(Eiter_2site))/L,abs(Eiter_2site(:)-E0_exact)+eps,'LineWidth',1);
set(gca,'XScale','Linear','YScale','log','FontSize',13,'LineWidth',1);
