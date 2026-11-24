clear; clc; 
% close all;

%% Step 1. loading data 
% ------>>>>>> read .nc file
file_mslp = 'mil0021_EH5_DJF_MSLP.nc'; 

% ------>>>>>> use 'ncdisp' to display all the data information in the .nc file
% you will see the names and dimensions of all the variables
% ncdisp(file_mslp); 

% ------>>>>>> use 'ncread' to read the variables you want to use
mslp_lat = ncread(file_mslp,'lat'); 
mslp_lon_orig = ncread(file_mslp,'lon');   
mslp_data = ncread(file_mslp,'slp'); %lon*lat*time
mslp_data = double(mslp_data)/100;   % convert Pa to hPa;
% Duration of model run
mslp_year = (800:800+length(mslp_data)-1)';    

%% Step 2. prepare the data for EOF analysis
% ------>>>>>> reorder the longitude from 0-360 to -180-180
mslp_lon = mod(mslp_lon_orig + 180, 360) - 180;
mslp_lon = [mslp_lon(49:end);mslp_lon(1:48)];
mslp_data_reorder = [mslp_data(49:end,:,:);mslp_data(1:48,:,:)];

% ------>>>>>> extract the North Atlantic area for EOF based on Hurrell(2003): 20N-70N, 90W-40E
lat_str = 20;  lat_end = 70;
lon_str = -90; lon_end = 40;
mslp_NAO = mslp_data_reorder(find(mslp_lon<=lon_str,1,'last'):find(mslp_lon>=lon_end,1,'first'),find(mslp_lat>=lat_end,1,'last'):find(mslp_lat<=lat_str,1,'first'),:); %lon*lat*time
mslp_NAO_lon = mslp_lon(find(mslp_lon<=lon_str,1,'last'):find(mslp_lon>=lon_end,1,'first'));
mslp_NAO_lat = mslp_lat(find(mslp_lat>=lat_end,1,'last'):find(mslp_lat<=lat_str,1,'first'));

%% Step 3. EOF analysis
mslp_NAO = permute(mslp_NAO,[2 1 3]); % Permute data (lat,lon,time)
no_of_PCs = 100; % number of PCs to be calculated

% ********************** run EOF *********************
[COEFF,PCs,VARexp,LAMBDA,d_LAMBDA] = Rean_PCA_anomalies(mslp_NAO_lat,mslp_NAO_lon,mslp_NAO,no_of_PCs);
% ****************************************************

% normalize PC time series
PCs = zscore(PCs,0,1);

% ----->>>>> plot the explained variance of the first 30 modes 
figure; hold on;
plot(1:30,VARexp(1:30));
xlabel('mode'); ylabel('explained variance');

%% Step 4. Plot the EOF modes onto maps

% p-value for masking out non-significant data
p_val = 0.01;

% Make custom colorbar that goes from blue to red
b = [ones(10,1);linspace(1,0,11)'];
g = [linspace(0,.9,10)';linspace(1,0,11)'];
r = [linspace(0,.9,10)';ones(11,1)];
cmap = [r g b];

% Set limits on colorbar
MaxData = 5; MinData = -MaxData;

% set the area of plot 
lat_str_plot = 20;    lat_end_plot = 90;
lon_str_plot = -120;  lon_end_plot = 60;

% extract the data with the plot area
mslp_plot = mslp_data_reorder(find(mslp_lon<=lon_str_plot,1,'last'):find(mslp_lon>=lon_end_plot,1,'first'),find(mslp_lat<=lat_end_plot,1,'first'):find(mslp_lat<=lat_str_plot,1,'first'),:); %lon*lat*time
mslp_lon_plot = mslp_lon(find(mslp_lon<=lon_str_plot,1,'last'):find(mslp_lon>=lon_end_plot,1,'first'));
mslp_lat_plot = mslp_lat(find(mslp_lat<=lat_end_plot,1,'first'):find(mslp_lat<=lat_str_plot,1,'first'));

% plot the first 3 leading EOF modes (spatial patterns)
figure; hold on;

for i_EOF = 1:3 
    
    subplot(1,3,i_EOF); hold on;
    title(['DJF MSLP: EOF ',num2str(i_EOF), ' (',num2str(roundn(VARexp(i_EOF)*100,-2)),'%)']);
   
    % ************************************************
    % ----->>>>> plot the regression of PC timeseries with climate field data 
    [regr_PC_mslp] = Rean_regr(PCs(:,i_EOF),permute(mslp_plot,[2 1 3])); % Permute data (lat,lon,time)
    EOF_plot = regr_PC_mslp.slope;
    % Mask out non-significant data
    EOF_plot(regr_PC_mslp.p_val > p_val) = NaN;
    % ************************************************
    
    % set the axis of map
    axesm('lambert','MapLatLimit',[lat_str_plot lat_end_plot],'MapLonLimit',[lon_str_plot lon_end_plot]);
    axis off; framem off; gridm on; mlabel on; plabel on;
    setm(gca,'PLinelocation',20,'PLabellocation',20,'fontsize',7);
    setm(gca,'MLinelocation',30,'MLabellocation',30,'fontsize',7);

    % plot data
    surfm(mslp_lat_plot, mslp_lon_plot, EOF_plot);
    
    % set colormap
    colormap(cmap);
    % Set limits on color bar
    set(gca, 'CLim', [MinData, MaxData]);
    % Add color bar
    h_cl = colorbar('location','southoutside');
    % move the position of the color bar for a nice figure
    ax = gca; ax_pos = ax.Position;
    h_cl.Position(2) = h_cl.Position(2)-0.2;
    ax.Position = ax_pos;
    
    % Add coastlines
    coastlines = load('coastlines');
    geoshow(coastlines.coastlat,coastlines.coastlon,'displaytype','line','color',[0 0 0])

end
% set the size of the figure
set(gcf,'unit','centimeters','position',[7,5,25,10]);
% export figure
file_name_out = 'EOF_NAO_mslp_djf';
print(gcf,file_name_out,'-r300','-dtiff');

%% Step 5. Compare the station-based and EOF-based NAO index

% ----->>>>> calculate the station-based NAO (from BERN03_Ex1_NAO_example.m)

% Reykjavik 
% Extract site data for NAO
site_lon = 360-21.94; % longitude in degrees East in the model lon_mat grid
site_lat = 64.15;
% Interpolate and extract timeseries
[mslp_rey_djf] = site_interp_echamT31(site_lon,site_lat,permute(mslp_data,[3 2 1])); % permute model data to (time, lat, lon)

% Ponta Delgada
site_lon = 360-25.68; % longitude in degrees East in the model lon_mat grid
site_lat = 37.74;
% Interpolate and extract timeseries
[mslp_ponta_djf] = site_interp_echamT31(site_lon,site_lat,permute(mslp_data,[3 2 1])); % permute model data to (time, lat, lon)

% Calculate NAO timeseries. zscore normalizes the data
NAO_station_based_djf = zscore(mslp_ponta_djf)-zscore(mslp_rey_djf);
% Make smoothed version of NAO timeseries
NAO_station_based_djf_smooth = smoothdata(NAO_station_based_djf,'movmean',31); % 31-year moving mean

% ----->>>>> Plot the station-based and EOF-based NAO timeserie
figure; hold on;

% Plot the station-based NAO index
plot(mslp_year,NAO_station_based_djf,'color',[1 1 1]*0.5) % grey color for raw data

% Plot the EOF-based NAO index
plot(mslp_year,PCs(:,1),'r');

% add the correlation coefficient between the station-based and EOF-based NAO
corr_NAO_station_EOF = corr(NAO_station_based_djf,PCs(:,1));
text(810,4,strcat('r=',num2str(roundn(corr_NAO_station_EOF,-2))),'Fontsize',11);

% add legend for the figure
legend('NAO:station-based','NAO:EOF-based','Orientation','horizontal','Location','NorthOutside');

% Define limits for plot
axis([790 2110 -4.5 4.5]);

% set the size of the figure
set(gcf,'unit','centimeters','position',[7,5,25,5]);
% export figure
file_name_out = 'EOF_vs_station_NAO_mslp_djf';
print(gcf,file_name_out,'-r300','-dtiff');