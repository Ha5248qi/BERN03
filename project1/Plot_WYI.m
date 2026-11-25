function [] = Plot_WYI(WYI)

years = (1:length(WYI)) + 799; 

figure;
plot(years, WYI, 'Color', [0.7 0.7 0.7]); 
hold on;

plot(years, movmean(WYI, 30), 'b', 'LineWidth', 1.5);
yline(0, 'k--'); 
xlabel('Year');
ylabel('WYI Anomaly (m/s)');
title('Indian Summer Monsoon Index (Webster-Yang)');
legend('Annual Anomaly', '30-year Running Mean');

end