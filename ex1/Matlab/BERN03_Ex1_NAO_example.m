%% load ECHAM5 data
mslp_mil0021_djf = netcdf_simple('../mil0021_EH5_DJF_MSLP.nc','slp');
temp2_mil0021_djf = netcdf_simple('../mil0021_EH5_DJF_TEMP2.nc','2t');
tprec_mil0021_djf = netcdf_simple('../mil0021_EH5_DJF_TPREC.nc','precip');
% Load lat-lon grid in matrix format
load lon_lat_mat.mat

%% NAO timeseries
% Extract site data for NAO
% Reykjavik 
site_lon = 360-21.94; % longitude in degrees East in the model lon_mat grid
site_lat = 64.15;
% Interpolate and extract timeseries
[mslp_rey_djf] = site_interp_echamT31(site_lon,site_lat,permute(mslp_mil0021_djf,[3 2 1])); % permute model data to (time, lat, lon)
%
% Ponta Delgada
site_lon = 360-25.68; % longitude in degrees East in the model lon_mat grid
site_lat = 37.74;
% Interpolate and extract timeseries
[mslp_ponta_djf] = site_interp_echamT31(site_lon,site_lat,permute(mslp_mil0021_djf,[3 2 1])); % permute model data to (time, lat, lon)
% Calculate NAO timeseries. zscore normalizes the data
NAO_station_based_djf = zscore(mslp_ponta_djf)-zscore(mslp_rey_djf);

%% Regression of NAO timeseries with DJF climate field data
% Regrssion with MSLP
[regr_NAO_mslp_djf] = Rean_regr(NAO_station_based_djf,permute(mslp_mil0021_djf,[2 1 3])); % Permute data (lat,lon,time)
% Regrssion with TEMP2
[regr_NAO_temp2_djf] = Rean_regr(NAO_station_based_djf,permute(temp2_mil0021_djf,[2 1 3])); % Permute data (lat,lon,time)
% Regression with TPREC
[regr_NAO_tprec_djf] = Rean_regr(NAO_station_based_djf,permute(tprec_mil0021_djf,[2 1 3])); % Permute data (lat,lon,time)

%% Plot NAO timeseries

% Define timescale
time_model = 800:2101; % Duration of model run

% Make smoothed version of NAO timeseries
NAO_station_based_djf_smooth = smoothdata(NAO_station_based_djf,'movmean',31); % 31-year moving mean

% Header to make plot look nice
PSize = [20,16];
h_fig = figure(1);
set(h_fig,'PaperUnits','centimeters','PaperSize',PSize,'PaperPosition',[0 0 PSize],'Color','w','InvertHardCopy','off');

% Make plot
% Use subplot to shape figure axis
subplot(3,3,[1 3])
% Plot raw timeseries
plot(time_model,NAO_station_based_djf,'color',[1 1 1]*0.5) % grey color for raw data
hold on % don't overwite next plot
% Plot smoothed data
plot(time_model,NAO_station_based_djf_smooth,'k','linewidth',1) % black color for raw data, with a thicker line

% Define limits for plot
axis([790 2110 -4.5 4.5])

% Add grid to plot
grid on

% Title and axis labels
title('mil0021 NAO DJF')
xlabel('Year')
ylabel('NAO index [norm.]')

% Print figure
print(gcf, '-depsc2', '-r300', 'mil0021_NAO_djf.eps');
close(h_fig)

%% Plot maps of correlation between climate field data and NAO timeseries

% p-value for masking out non-significant correlations
p_val = 0.01;

% Closing the gap in the plot
lon_mat(:,end)=lon_mat(:,1);

% Make custom colorbar that goes from blue to red
b = [ones(10,1);linspace(1,0,11)'];
g = [linspace(0,.9,10)';linspace(1,0,11)'];
r = [linspace(0,.9,10)';ones(11,1)];
cmap = [r g b];
colormap(cmap);

% Set limets on colorbar
MaxData = 1;
MinData = -MaxData;

% Header to make plot look nice
PSize = [20,16];
h_fig = figure(1);
set(h_fig,'PaperUnits','centimeters','PaperSize',PSize,'PaperPosition',[0 0 PSize],'Color','w','InvertHardCopy','off');

% ******************************* NAO corr MSLP **************************
% Define size of subplot
subplot('position', [0.02 0.5 0.3 0.4]); % [left bottom width height]
% Use Lambert projection and show North Atlantic region
axesm('lambert','MapLatLimit',[20 90],'MapLonLimit',[-130 70])
% Options for map
axis off; framem off; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-90,'MLabelLocation',[-150 -90 -30 0 30 90 150],'fontsize',7)

% Prepare data for plotting
plot_data = double(squeeze(regr_NAO_mslp_djf.cor));
% Mask out non-significant correlations
plot_data(regr_NAO_mslp_djf.p_val > p_val) = NaN;
% plot data
surfm(lat_mat, lon_mat, plot_data)

% Set limets on color bar
set(gca, 'CLim', [MinData, MaxData]);
% Add color bar
colorbar('location','southoutside')

% Add coastlines
coastlines = load('coastlines');
geoshow(coastlines.coastlat,coastlines.coastlon,'displaytype','line','color',[0 0 0])
% Title
title('NAO corr. MSLP DJF')

% ******************************* NAO corr TEMP2 **************************
% Define size of subplot
subplot('position', [0.34 0.5 0.3 0.4]); % [left bottom width height]
% Use Lambert projection and show North Atlantic region
axesm('lambert','MapLatLimit',[20 90],'MapLonLimit',[-130 70])
% Options for map
axis off; framem off; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-90,'MLabelLocation',[-150 -90 -30 0 30 90 150],'fontsize',7)

% Prepare data for plotting
plot_data = double(squeeze(regr_NAO_temp2_djf.cor));
% Mask out non-significant correlations
plot_data(regr_NAO_temp2_djf.p_val > p_val) = NaN;
% plot data
surfm(lat_mat, lon_mat, plot_data)

% Set limets on color bar
set(gca, 'CLim', [MinData, MaxData]);
% Add color bar
colorbar('location','southoutside')

% Add coastlines
coastlines = load('coastlines');
geoshow(coastlines.coastlat,coastlines.coastlon,'displaytype','line','color',[0 0 0])
% Title
title('NAO corr. TEMP2 DJF')

% ******************************* NAO corr TPREC **************************
% Define size of subplot
subplot('position', [0.66 0.5 0.3 0.4]); % [left bottom width height]
% Use Lambert projection and show North Atlantic region
axesm('lambert','MapLatLimit',[20 90],'MapLonLimit',[-130 70])
% Options for map
axis off; framem off; gridm on; mlabel on; plabel on;
setm(gca,'MLabelParallel',-90,'MLabelLocation',[-150 -90 -30 0 30 90 150],'fontsize',7)

% Prepare data for plotting
plot_data = double(squeeze(regr_NAO_tprec_djf.cor));
% Mask out non-significant correlations
plot_data(regr_NAO_tprec_djf.p_val > p_val) = NaN;
% plot data
surfm(lat_mat, lon_mat, plot_data)

% Set limets on color bar
set(gca, 'CLim', [MinData, MaxData]);
% Add color bar
colorbar('location','southoutside')

% Add coastlines
coastlines = load('coastlines');
geoshow(coastlines.coastlat,coastlines.coastlon,'displaytype','line','color',[0 0 0])
% Title
title('NAO corr. TPREC DJF')

% Print figure
print(gcf, '-depsc2', '-r300', 'mil0021_NAO_corr_MSLP_TEMP2_TPREC_djf.eps');
close(h_fig)
