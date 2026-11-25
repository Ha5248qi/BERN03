
%% Data extraction


u850_mil0021_JJA = netcdf_simple('data/mil0021_JJA_DJF/mil0021_EH5_JJA_U850.nc','u');
u200_mil0021_JJA = netcdf_simple('data/mil0021_JJA_DJF/mil0021_EH5_JJA_U200.nc','u');
u850_mil0025_JJA = netcdf_simple('data/mil0025_JJA_DJF/mil0025_EH5_JJA_U850.nc','u');
u200_mil0025_JJA = netcdf_simple('data/mil0025_JJA_DJF/mil0025_EH5_JJA_U200.nc','u');
u850_mil0026_JJA = netcdf_simple('data/mil0026_JJA_DJF/mil0026_EH5_JJA_U850.nc','u');
u200_mil0026_JJA = netcdf_simple('data/mil0026_JJA_DJF/mil0026_EH5_JJA_U200.nc','u');


lon = netcdf_simple('data/mil0021_JJA_DJF/mil0021_EH5_JJA_U850.nc','lon');
lat = netcdf_simple('data/mil0021_JJA_DJF/mil0021_EH5_JJA_U850.nc','lat');

%% Webster index

WYI_0021_anom = webster_index(u850_mil0021_JJA,u200_mil0021_JJA)
WYI_0025_anom = webster_index(u850_mil0025_JJA,u200_mil0025_JJA)
WYI_0026_anom = webster_index(u850_mil0026_JJA,u200_mil0026_JJA)



%% Plot separately

 Plot_WYI(WYI_0021_anom)

 Plot_WYI(WYI_0025_anom)

 Plot_WYI(WYI_0026_anom)

 %% comparison plot
 years = (1:length(WYI_0021_anom)) + 799; 

 figure;
 hold on;

 plot(years, movmean(WYI_0021_anom, 30), 'b', 'LineWidth', 1.5);
 plot(years, movmean(WYI_0025_anom, 30), 'r', 'LineWidth', 1.5);
 plot(years, movmean(WYI_0026_anom, 30), 'Color', [0 0.6 0], 'LineWidth', 1.5);

 yline(0, 'k--'); 

 xlabel('Year');
 ylabel('WYI Anomaly (m/s)');
 title('Comparison of Indian Summer Monsoon Index (30-year Running Mean)');


 legend('mil0021', 'mil0025', 'mil0026 ', 'Location', 'best');
 grid on;
 hold off;

 %% Variance

 var_annual = [var(WYI_0021_anom), var(WYI_0025_anom), var(WYI_0026_anom)];
 var_30yr = [var(movmean(WYI_0021_anom, 30)), var(movmean(WYI_0025_anom, 30)), var(movmean(WYI_0026_anom, 30))];

 
 fprintf('0021: Annual Var = %.3f , 30-yr Var = %.4f\n', var_annual(1),  var_30yr(1));
 fprintf('0025: Annual Var = %.3f , 30-yr Var = %.4f\n', var_annual(2),  var_30yr(2));
 fprintf('0026: Annual Var = %.3f , 30-yr Var = %.4f\n', var_annual(3), var_30yr(3));



 %% Volcanic affection

 volcdata = readmatrix('data/volc_data','Filetype','text','CommentStyle','#');

 % Volcanic aerosols
 
 Volc_annual(:,1) = (790:2009)';
 Volc_annual(:,2) = mean(reshape(volcdata(1:end-1,2),36,[]),1)';  % 30°N–90°N
 Volc_annual(:,3) = mean(reshape(volcdata(1:end-1,4),36,[]),1)';  % 0°–30°N
 Volc_annual(:,4) = mean(reshape(volcdata(1:end-1,6),36,[]),1)';  % 0°–30°S
 Volc_annual(:,5) = mean(reshape(volcdata(1:end-1,8),36,[]),1)';  % 30°S–90°S

 % Plot the four latitude bands in separate subplots (North to South)
 % figure;
 % 
 % subplot(4,1,1);
 % plot(Volc_annual(:,1), Volc_annual(:,2), 'b');
 % axis([800 2100 0 0.7]);
 % title('30°N–90°N');
 % ylabel('AOD');
 % grid on;
 % 
 % subplot(4,1,2);
 % plot(Volc_annual(:,1), Volc_annual(:,3), 'r');
 % axis([800 2100 0 0.7]);
 % title('0°–30°N');
 % ylabel('AOD');
 % grid on;
 % 
 % subplot(4,1,3);
 % plot(Volc_annual(:,1), Volc_annual(:,4), 'm');
 % axis([800 2100 0 0.7]);
 % title('0°–30°S');
 % ylabel('AOD');
 % grid on;
 % 
 % subplot(4,1,4);
 % plot(Volc_annual(:,1), Volc_annual(:,5), 'g');
 % axis([800 2100 0 0.7]);
 % title('30°S–90°S');
 % xlabel('Year');
 % ylabel('AOD');
 % grid on;


 %% Eruption Identification
 n_years = 1220;
 n_per_year = 36;

 AOD_0_30N = volcdata(1:n_years*n_per_year, 4);
 AOD_0_30S = volcdata(1:n_years*n_per_year, 6);

 % Reshape into [36 time steps, 1220 years]
 AOD_0_30N = reshape(AOD_0_30N, n_per_year, n_years);
 AOD_0_30S = reshape(AOD_0_30S, n_per_year, n_years);


 AOD_0_30N_annual = mean(AOD_0_30N, 1)';
 AOD_0_30S_annual = mean(AOD_0_30S, 1)';
 years_volc = (790:2009)';

 AOD_0_30N_anom = AOD_0_30N_annual - mean(AOD_0_30N_annual);
 AOD_0_30S_anom = AOD_0_30S_annual - mean(AOD_0_30S_annual);

 %% 4. Plot annual AOD anomalies

 % figure;
 % 
 % subplot(2,1,1);
 % plot(years_volc, AOD_0_30N_anom, 'b');
 % yline(0, '--k');
 % xlabel('Year');
 % ylabel('AOD anomaly');
 % title('AOD Anomaly (0-30N)');
 % grid on;
 % xlim([790 2010]);
 % 
 % subplot(2,1,2);
 % plot(years_volc, AOD_0_30S_anom, 'r');
 % yline(0, '--k');
 % xlabel('Year');
 % ylabel('AOD anomaly');
 % title('AOD Anomaly (0-30S)');
 % grid on;
 % xlim([790 2010]);
 % 
 % sgtitle('Tropical and Subtropical Aerosol Optical Depth Anomalies');

 %% combined plot
 % 
 % figure;
 % 
 % plot(years_volc, AOD_0_30N_anom, 'b', 'LineWidth', 1);
 % hold on;
 % plot(years_volc, AOD_0_30S_anom, 'r', 'LineWidth', 1);
 % 
 % yline(0, '--k');
 % 
 % xlabel('Year');
 % ylabel('AOD anomaly');
 % title('Tropical and Subtropical AOD Anomalies');
 % legend('0-30N', '0-30S', 'Location', 'best');
 % 
 % grid on;
 % xlim([790 2010]);
 % hold off;

 %% Select top eruption years

 % Restrict eruption search range so WYI covers [-10, +5] years around each event
 ok = (years_volc >= min(years) + 10) & (years_volc <= max(years) - 5);
 yr = years_volc(ok);

 % Tropical eruptions (mean of 0-30N and 0-30S)
 AOD_trop = 0.5 * (AOD_0_30N_annual(ok) + AOD_0_30S_annual(ok));
 [~, erup_trop] = findpeaks(AOD_trop, yr, 'MinPeakDistance', 10, 'SortStr', 'descend');
 erup_trop = erup_trop(1:8); % Top 8

 AOD_30_90N = Volc_annual(ok, 2); % 30-90N (target band)
 AOD_0_30S  = Volc_annual(ok, 4); % 0-30S (to filter out tropical eruptions)
 [~, locs]  = findpeaks(AOD_30_90N, 'MinPeakDistance', 10, 'MinPeakHeight', 0.02, 'SortStr', 'descend');
 locs       = locs(AOD_0_30S(locs) < 0.02);
 erup_nh    = yr(locs(1:min(8, end)));

 disp('Tropical eruption years:');    disp(erup_trop');
 disp('NH high-lat eruption years:'); disp(erup_nh');


 %% Superposed Epoch Analysis (SEA) for WYI

 lag = (-10:5)'; % 10 years before to 5 years after eruption
 WYI_all = [WYI_0021_anom(:), WYI_0025_anom(:), WYI_0026_anom(:)]; % [time x 3 members]

 %Tropical eruptions SEA
 sea_trop = [];
 for i = 1:length(erup_trop)
     idx  = find(years == erup_trop(i));
     win  = WYI_all(idx + lag, :);
     base = mean(WYI_all(idx + (-10:-1), :), 1); % Pre-eruption 10-yr mean
     sea_trop = [sea_trop, win - base];
 end

 % NH high-latitude eruptions SEA
 sea_nh = [];
 for i = 1:length(erup_nh)
     idx  = find(years == erup_nh(i));
     win  = WYI_all(idx + lag, :);
     base = mean(WYI_all(idx + (-10:-1), :), 1); % Pre-eruption 10-yr mean
     sea_nh = [sea_nh, win - base];
 end

 % Mean and standard error (SE) across events and ensemble members
 mean_trop = mean(sea_trop, 2);
 se_trop   = std(sea_trop, 0, 2) / sqrt(size(sea_trop, 2));

 mean_nh = mean(sea_nh, 2);
 se_nh   = std(sea_nh, 0, 2) / sqrt(size(sea_nh, 2));

 % Plot SEA comparison
 figure;
 hold on;
 errorbar(lag, mean_trop, se_trop, 'b-o', 'LineWidth', 1.5);
 errorbar(lag, mean_nh,   se_nh,   'r-s', 'LineWidth', 1.5);
 yline(0, 'k--');
 xline(0, 'k:', 'Eruption Year', 'LabelVerticalAlignment', 'bottom');

 xlabel('Lag (Years relative to eruption)');
 ylabel('WYI Anomaly (m/s)');
 title('SEA of Indian Summer Monsoon (WYI) to Volcanic Eruptions');
 legend('Tropical Eruptions', 'NH High-Lat Eruptions', 'Location', 'best');
 grid on;
 hold off;



 %% Monte Carlo Significance Test (95% CI)

 n_mc = 2000;
 valid_idx = find(years >= min(years)+10 & years <= max(years)-5);

 figure('Color', 'w', 'Position', [100, 100, 920, 400]);
 means  = {mean_trop, mean_nh};
 n_ev   = [length(erup_trop), length(erup_nh)];
 titles = {'Tropical Eruptions', 'NH High-Lat Eruptions'};

 for g = 1:2
     mc_means = zeros(length(lag), n_mc);
     for k = 1:n_mc
         r_idx = valid_idx(randperm(length(valid_idx), n_ev(g)));
         tmp = [];
         for i = 1:n_ev(g)
             win  = WYI_all(r_idx(i) + lag, :);
             base = mean(WYI_all(r_idx(i) + (-10:-1), :), 1);
             tmp  = [tmp, win - base];
         end
         mc_means(:, k) = mean(tmp, 2);
     end

     %ci_low  = prctile(mc_means, 2.5, 2);
     %ci_high = prctile(mc_means, 97.5, 2);   % two trail
     ci_low = prctile(mc_means, 5, 2); % one trail

     subplot(1, 2, g); hold on;
     plot(lag, means{g}, 'b-o', 'LineWidth', 1.5); 
     plot(lag, ci_low,  'k--', 'LineWidth', 1.1); 
     %plot(lag, ci_high, 'k--', 'LineWidth', 1.1); 
     yline(0, 'k:'); xline(0, 'k:');
     xlim([-10, 5]);
     title(titles{g});
     xlabel('Lag (Years)');
     ylabel('WYI Anomaly (m/s)');
     grid on;
 end
