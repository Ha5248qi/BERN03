function [data] = netcdf_simple(file,varname)
% Open file.
ncid = netcdf.open(file,'NC_NOWRITE');
% Get variable ID of the variable, given its name.
varid = netcdf.inqVarID(ncid,varname);
% Get the value of the variable, given its ID.
data_read = netcdf.getVar(ncid,varid,'double');
data = squeeze(data_read);