errors = [0.0675,0.0329,0.0192,0.0124,0.0087,0.0064];
Nstep = [20,30,40,50,60,70];
deltaT = 50* ones(1, 6)./Nstep;
deltaT_2 = deltaT .* deltaT; 
deltaT_3 = deltaT_2.*deltaT;


p = polyfit(log(deltaT), log(errors), 1)
y1 = polyval(p, log(deltaT));

figure;
hold on 
plot(deltaT,errors,'o','LineWidth',1);
plot(deltaT, exp(y1), '--','LineWidth',1);

set(gca, 'Xscale', 'log', 'Yscale', 'log', 'FontSize',13,'LineWidth',1)
grid on;
title('Purification error');
xlabel('\Delta t')
ylabel('Maximum relative error');
