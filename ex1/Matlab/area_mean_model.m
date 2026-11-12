function [outdata] = area_mean_model(data,lat_mat,lon_mat,lon_east,lon_west,lat_north,lat_south)
% function calculates area weighted mean given lon and lat limits and grid
% data has format (time,lat,lon)
t_length = size(data,1);
lat_length  = size(data,2);
lon_length  = size(data,3);
% for area weighting
e=0.08181921804834; % Earth's eccentricity
Atot=510072000; % total area [km2]
a=6378.1370; % semi-major axis [km]
% Formula for area weighting A = (a^2*(1-e^2)*sin(alpha)^2*cos(phi))/(1-e^2*sin(phi)^2)^2 
phi = lat_mat*pi/180;
alpha = 180/lat_length;
A = (a^2*(1-e^2)*sin(alpha)^2*cos(phi))./((1-e^2*sin(phi).*sin(phi)).*(1-e^2*sin(phi).*sin(phi)));
% area_w = A/(sum(sum(A)))*lat_length*lon_length; % area weighting matrix
area_w = sqrt(cos(phi)); % simple area weighting
area_w_mat = repmat(reshape(area_w,1,lat_length,lon_length),[t_length,1,1]);
% cut out area
data_out = nan(size(data));
%data_out(:,lat_mat>=lat_south & lat_mat<=lat_north & lon_mat<=lon_east | lat_mat>=lat_south & lat_mat<=lat_north & lon_mat>=lon_west) = data ... 
%    (1:t_length,lat_mat>=lat_south & lat_mat<=lat_north & lon_mat<=lon_east | lat_mat>=lat_south & lat_mat<=lat_north &lon_mat>=lon_west);
data_out(:,lat_mat>=lat_south & lat_mat<=lat_north & lon_mat<=lon_east & lat_mat>=lat_south & lat_mat<=lat_north & lon_mat>=lon_west) = data ... 
   (1:t_length,lat_mat>=lat_south & lat_mat<=lat_north & lon_mat<=lon_east & lat_mat>=lat_south & lat_mat<=lat_north &lon_mat>=lon_west);
%data(:,lon_mat >= lon_limE) = NaN;
%data(:,lon_mat <= lon_limW) = NaN;
%data(:,lat_mat >= lat_limN) = NaN;
%data(:,lat_mat <= lat_limS) = NaN;
% mean values weighted by area
data_out = squeeze(nanmean(nanmean(data_out.*area_w_mat,3),2));
outdata = data_out;