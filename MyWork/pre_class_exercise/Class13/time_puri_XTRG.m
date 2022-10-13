
%Elapsed time: 10.3s, CPU time: 385.2s, Avg # of cores: 37.4
%Elapsed time: 17.25s, CPU time: 674.9s, Avg # of cores: 39.12
%Elapsed time: 24.51s, CPU time: 970.3s, Avg # of cores: 39.59

% Elapsed time: 4.69s, CPU time: 181s, Avg # of cores: 38.58
% Elapsed time: 8.95s, CPU time: 351.5s, Avg # of cores: 39.27
% Elapsed time: 13.11s, CPU time: 514.9s, Avg # of cores: 39.28

puri = [17.14,27.15,39.19,46.57,57.04];
XTRG =[8.981,16.9,24.75,33.1,41.22];
Nsweep = (1:5);
figure; 
hold on 
plot(Nsweep, puri, '-x','LineWidth',1);
plot(Nsweep, XTRG, '-x','LineWidth',1);

hold off
set(gca,'XScale','log','YScale','log','FontSize',13,'LineWidth',1)
grid on;
%title('Purification');
xlabel('Nsweep')
ylabel('Elapsed time [s]');