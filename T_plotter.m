figure(1)
plot(x,423*ones(size(x))); grid on; hold on;
xlim([0, max(x)]);
xlabel('Length [m]');
ylabel('Temperature [K]');
title('Transient behavior of mean wall temperature')



NPLOT = 4;
tPLOT = [25 50 100];

for i = 1:length(TMID(1,:))
    if max(i*Dt == tPLOT)
        T_plot = TMID(:,i)';
        plot(x,T_plot);
        tPLOT = [tPLOT, time];
    end
end
text(0.005, 425, 't = 0', 'Color', [0, 0.4470, 0.7410])
text(0.005, 433, 't = 25', 'Color', [0.8500, 0.3250, 0.0980])
text(0.005, 438.5, 't = 50', 'Color',[0.9290, 0.6940, 0.1250])
text(0.04, 445, 't = 100', 'Color',[0.4940, 0.1840, 0.5560])
