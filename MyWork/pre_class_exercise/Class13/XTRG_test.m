clear

% system parameter
J = -1; % coupling strength
L = 20; % number of sites in a chain

% imaginary-time evolution parameters
beta = 100; % inverse temperature
Nstep = 40; % # of iterations
Nkeep = 40; % maximum bond dimension
Nsweep = 2; % # of round-trip sweeps in the variational multiplication

% Local operators
[S,I] = getLocalSpace('Spin',1/2);

Hloc = cell(4,4);
Hloc(:) = {zeros(size(I))};
Hloc{1,1} = I;
Hloc{2,1} = S(:,:,1);
Hloc{3,1} = S(:,:,3);
Hloc{4,2} = J*S(:,:,1)';
Hloc{4,3} = J*S(:,:,3)';
Hloc{end,end} = I;
Hloc = cell2mat(reshape(Hloc,[1 1 size(Hloc,1) size(Hloc,2)]));

% full chain
Hs = cell(1,L);
Hs(:) = {Hloc};
Hs{1} = Hs{1}(:,:,end,:); % choose the last components of the left leg
Hs{end} = Hs{end}(:,:,:,1); % choose the first components of the right leg

tmax = beta;
dt = tmax/(2^Nstep);
for Nsweep =(1:5)
    [taus,lnZs,rho] = XTRG_Ex(Hs,dt,tmax,Nkeep,Nsweep);
end
% figure;
% lgd = cell(1, 6);
% step = 1;
% hold on;
% for Nkeep = (10:10:60)
%     [taus,lnZs,rho] = XTRG_Ex(Hs,dt,tmax,Nkeep,Nsweep);
%     betas = taus;
% % exact values
%     epsk = J*cos((1:L).'*(pi/(L+1)));
%     lnZexact = sum(log(1+exp(epsk.*(-betas))),1);
% 
%     plot(betas,abs(lnZs./lnZexact-1),'-x','LineWidth',1);
%     err = abs(lnZs./lnZexact-1);
%     max(err)
%     lgd{step} = strcat('Nkeep=', int2str(Nkeep));
%     step  = step+1;
% end 
% 
% grid on;
% set(gca,'XScale','log','YScale','log','FontSize',13,'LineWidth',1)
% 
% %axis([10^-11, 10^2, 10^-16, 10^-2]);
% title('XTRG');
% xlabel('\beta')
% ylabel('Relative error in the partition function');
% legend(lgd, 'Location','southeast');
% AB first
% Elapsed time: 10.19s, CPU time: 317.5s, Avg # of cores: 31.14
% Elapsed time: 18.21s, CPU time: 568.7s, Avg # of cores: 31.23
% Elapsed time: 25.91s, CPU time: 827.8s, Avg # of cores: 31.95

% Elapsed time: 4.69s, CPU time: 181s, Avg # of cores: 38.58
% Elapsed time: 8.95s, CPU time: 351.5s, Avg # of cores: 39.27
% Elapsed time: 13.11s, CPU time: 514.9s, Avg # of cores: 39.28