clear

% system parameter
J = -1; % coupling strength
L = 20; % number of sites in a chain

% imaginary-time evolution parameters
beta = 100; % inverse temperature
tmax = beta/2;
Nstep = 40; % # of iterations
%dt = tmax/Nstep; % time step
Nkeep = 40; % maximum bond dimension
Nsweep = 2; % # of round-trip sweeps in the variational multiplication

% Local operators
[S,I] = getLocalSpace('Spin',1/2);
% nearest-neighbor interaction terms
Hs = cell(1,L-1);
Hs(:) = {J*contract(S(:,:,[1 3]),3,3, ...
    permute(conj(S(:,:,[1 3])),[2 1 3]),3,3)};
dt = tmax/Nstep; % time step
for Nsweep =(1:5)
    [taus,lnZs,MPO] = finT_puri_Ex (Hs,dt,tmax,Nkeep,Nsweep);
end
% betas = taus*2;
% 
% % exact values
% epsk = J*cos((1:L).'*(pi/(L+1)));
% lnZexact = sum(log(1+exp(epsk.*(-betas))),1);
% 
% %figure;
%     plot(betas, lnZs./lnZexact-1,'-x','LineWidth',1);
% figure;
% lgd = cell(1, 6);
% step = 1;
% hold on;
% for Nstep = (20:10:70) % # of iterations
%     dt = tmax/Nstep; % time step
%     [taus,lnZs,MPO] = finT_puri_Ex (Hs,dt,tmax,Nkeep,Nsweep);
% 
%     betas = taus*2;
%     
%     % exact values
%     epsk = J*cos((1:L).'*(pi/(L+1)));
%     lnZexact = sum(log(1+exp(epsk.*(-betas))),1);
% 
% %figure;
%     plot(betas, lnZs./lnZexact-1,'-x','LineWidth',1);
%     err = max(lnZs./lnZexact-1);
%     err(end)
%     lgd{step} = strcat('Nstep=', int2str(Nstep));
%     step  = step+1;
% end
% hold off;
% set(gca,'FontSize',13,'LineWidth',1)
% grid on;
% title('Purification');
% xlabel('\beta')
% ylabel('Relative error in the partition function');
% legend(lgd, 'Location','southeast')
% Contract AB first
% Elapsed time: 10.93s, CPU time: 413.3s, Avg # of cores: 37.8
% Elapsed time: 17.77s, CPU time: 698s, Avg # of cores: 39.27
% Elapsed time: 24.73s, CPU time: 1004s, Avg # of cores: 40.62

% NO
%Elapsed time: 10.3s, CPU time: 385.2s, Avg # of cores: 37.4
%Elapsed time: 17.25s, CPU time: 674.9s, Avg # of cores: 39.12
%Elapsed time: 24.51s, CPU time: 970.3s, Avg # of cores: 39.59
