%% load ECHAM5 data
% mil0021
mslp_mil0021_djf = netcdf_simple('mil0021_EH5_DJF_MSLP.nc','slp');
temp2_mil0021_djf = netcdf_simple('mil0021_EH5_DJF_TEMP2.nc','2t');
tprec_mil0021_djf = netcdf_simple('mil0021_EH5_DJF_TPREC.nc','precip');
% mil0025
mslp_mil0025_djf = netcdf_simple('mil0025_EH5_DJF_MSLP.nc','slp');
temp2_mil0025_djf = netcdf_simple('mil0025_EH5_DJF_TEMP2.nc','2t');
tprec_mil0025_djf = netcdf_simple('mil0025_EH5_DJF_TPREC.nc','precip');
% mil0026
mslp_mil0026_djf = netcdf_simple('mil0026_EH5_DJF_MSLP.nc','slp');
temp2_mil0026_djf = netcdf_simple('mil0026_EH5_DJF_TEMP2.nc','2t');
tprec_mil0026_djf = netcdf_simple('mil0026_EH5_DJF_TPREC.nc','precip');
% Load lat-lon grid in matrix format
load lon_lat_mat.mat

%%
% Annual forcing data
% TSI 790-2101
TSI_annual(:,1) = mean(reshape(bard25tsi2000plus(:,2),12,[]),1)';
TSI_annual(:,2) = mean(reshape(bard25tsi2000plus(:,3),12,[]),1)';

plot(TSI_annual(:,1),TSI_annual(:,2),'r')
axis([800 2100 1363 1370])
%%
% Volcanic aerosols
Volc_annual(:,1) = (790:2009)';
Volc_annual(:,2) = mean(reshape(volcdata(1:end-1,4),36,[]),1)';
Volc_annual(:,3) = mean(reshape(volcdata(1:end-1,7),36,[]),1)';
Volc_annual(:,4) = mean(reshape(volcdata(1:end-1,10),36,[]),1)';
Volc_annual(:,5) = mean(reshape(volcdata(1:end-1,13),36,[]),1)';

plot(Volc_annual(:,1),Volc_annual(:,2),'b');
hold on
plot(Volc_annual(:,1),Volc_annual(:,3),'r');
plot(Volc_annual(:,1),Volc_annual(:,4),'m');
plot(Volc_annual(:,1),Volc_annual(:,5),'g');
axis([800 2100 0 0.7])

%% Plot forcings
% CO2
subplot(3,3,[1 3])
plot(E2CO28002100(:,1),E2CO28002100(:,2),'k')
hold on
plot(E2CO28002100(:,1),E2CO28002100(:,3),'k')
plot(E2CO28002100(:,1),E2CO28002100(:,4),'k')
grid on
axis([800 2100 270 290])
title('CO2 [PPM]')

% Volc
subplot(3,3,[4 6])
plot(Volc_annual(:,1),Volc_annual(:,2),'b');
hold on
plot(Volc_annual(:,1),Volc_annual(:,3),'r');
plot(Volc_annual(:,1),Volc_annual(:,4),'m');
plot(Volc_annual(:,1),Volc_annual(:,5),'g');
grid on
axis([800 2100 0 0.7])
title('Volcanic aerosols [AOD]')

% TSI
subplot(3,3,[7 9])
plot(TSI_annual(:,1),TSI_annual(:,2),'r')
grid on
axis([800 2100 1363 1370])
title('Total solar irradiance [W/m2]')

xlabel('Year')

%%
% NAO
% NAO timeseries
% Extract site data for NAO
% Reykjavik 
site_lon = 360-21.94; % longitude in degrees East in the model lon_mat grid
site_lat = 64.15;
% Interpolate and extract timeseries
[mslp_mil0021_rey_djf] = site_interp_echamT31(site_lon,site_lat,permute(mslp_mil0021_djf,[3 2 1])); % permute model data to (time, lat, lon)
[mslp_mil0025_rey_djf] = site_interp_echamT31(site_lon,site_lat,permute(mslp_mil0025_djf,[3 2 1])); % permute model data to (time, lat, lon)
[mslp_mil0026_rey_djf] = site_interp_echamT31(site_lon,site_lat,permute(mslp_mil0026_djf,[3 2 1])); % permute model data to (time, lat, lon)
%
% Ponta Delgada
site_lon = 360-25.68; % longitude in degrees East in the model lon_mat grid
site_lat = 37.74;
% Interpolate and extract timeseries
[mslp_ponta_mil0021_djf] = site_interp_echamT31(site_lon,site_lat,permute(mslp_mil0021_djf,[3 2 1])); % permute model data to (time, lat, lon)
[mslp_ponta_mil0025_djf] = site_interp_echamT31(site_lon,site_lat,permute(mslp_mil0025_djf,[3 2 1])); % permute model data to (time, lat, lon)
[mslp_ponta_mil0026_djf] = site_interp_echamT31(site_lon,site_lat,permute(mslp_mil0026_djf,[3 2 1])); % permute model data to (time, lat, lon)

% Calculate NAO timeseries. zscore normalizes the data
NAO_station_based_mil0021_djf = zscore(mslp_ponta_mil0021_djf)-zscore(mslp_mil0021_rey_djf);
NAO_station_based_mil0025_djf = zscore(mslp_ponta_mil0025_djf)-zscore(mslp_mil0025_rey_djf);
NAO_station_based_mil0026_djf = zscore(mslp_ponta_mil0026_djf)-zscore(mslp_mil0026_rey_djf);

%%
% Define timescale
time_model = 800:2101; % Duration of model run

% Make smoothed version of NAO timeseries
NAO_station_based_mil0021_djf_smooth = smoothdata(NAO_station_based_mil0021_djf,'movmean',31); % 31-year moving mean
NAO_station_based_mil0025_djf_smooth = smoothdata(NAO_station_based_mil0025_djf,'movmean',31); % 31-year moving mean
NAO_station_based_mil0026_djf_smooth = smoothdata(NAO_station_based_mil0026_djf,'movmean',31); % 31-year moving mean

% Header to make plot look nice
PSize = [20,16];
h_fig = figure(2);
set(h_fig,'PaperUnits','centimeters','PaperSize',PSize,'PaperPosition',[0 0 PSize],'Color','w','InvertHardCopy','off');

% Make plot
% Use subplot to shape figure axis
subplot(3,3,[1 3])
% Plot raw timeseries
plot(time_model,NAO_station_based_mil0021_djf,'color',[1 1 1]*0.4) % grey color for raw data
hold on % don't overwite next plot
plot(time_model,NAO_station_based_mil0025_djf,'color',[1 1 1]*0.6) % grey color for raw data
plot(time_model,NAO_station_based_mil0026_djf,'color',[1 1 1]*0.8) % grey color for raw data
% Define limits for plot
axis([790 2110 -4.5 4.5])
grid on
ylabel('NAO index [norm.]')
title('E2 NAO DJF')

% Use subplot to shape figure axis
subplot(3,3,[4 6])
% Plot smoothed data
plot(time_model,NAO_station_based_mil0021_djf_smooth,'color',[1 1 1]*0.4,'linewidth',1) %  with a thicker line
hold on % don't overwite next plot
plot(time_model,NAO_station_based_mil0025_djf_smooth,'color',[1 1 1]*0.6,'linewidth',1) %  with a thicker line
plot(time_model,NAO_station_based_mil0026_djf_smooth,'color',[1 1 1]*0.8,'linewidth',1) %  with a thicker line

% Define limits for plot
axis([790 2110 -2 1.5])
legend('mil0021','mil0025','mil0026','location','southeast')

% Add grid to plot
grid on

% Title and axis labels
xlabel('Year')
ylabel('NAO index [norm.]')

% Print figure
print(gcf, '-dtiff', '-r300', 'E2_NAO_djf.tiff');
close(h_fig)
