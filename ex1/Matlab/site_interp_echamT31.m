function [site_data] = site_interp_echamT31(site_lon,site_lat,data)
% The function site_interp extracts a time series of REMOiso output data
% in rotated coordinates loaded in the variable 'data'. The variable 'data'
% is assumed to have the format of (lon,lat,time).
% 'site_lon' and 'site_lat' is the geographcal coordinates for the site.
% The output is produced by interpolating the REMOiso data in the rotated
% grid to the site in question.
% site_lon = -21.95; % http://www.fallingrain.com/world/IC/10/Reykjavik.html
% site_lat = 64.15;
% lat, lon
%[lat] = netcdf_global('62041t42',2,'circ');
%[lon] = netcdf_global('62041t42',1,'circ');
%[lon_mat,lat_mat] = meshgrid(lon(1:128),lat(1:64));

load lon_lat_mat

% ZI = interp2(X,Y,Z,XI,YI,method)
for i = 1:size(data,1)
    site_data(i,1) = interp2(lon_mat,lat_mat,squeeze(data(i,:,:)),site_lon,site_lat,'*linear');
end