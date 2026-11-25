function [WYI_anom] = webster_index(u850,u200)

%% Research  area

lon = netcdf_simple('data/mil0021_JJA_DJF/mil0021_EH5_JJA_U850.nc','lon');
lat = netcdf_simple('data/mil0021_JJA_DJF/mil0021_EH5_JJA_U850.nc','lat');

%% webster index

lon_india=find(lon >=40 & lon <=110);
lat_india=find(lat >=0 & lat <=20);

u_shear = u850(lon_india, lat_india, :) - u200(lon_india, lat_india, :);
WYI = squeeze(mean(mean(u_shear, 1), 2));
WYI_anom = WYI - mean(WYI);

end