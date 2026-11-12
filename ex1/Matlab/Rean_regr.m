%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%                          Rean_regress.m                     %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% this function computes a single variable least squares regression on a field (e.g.Reananlysis)
% as a matrix operation, which is faster than running a loop with regress.m for all grid cells
%
% significance levels are computed from the F-Statistic
%
% !!!!     note that this function will ignore NaNs   !!!!
%
% input parameters
%       x: predictor (2-D vector, time series data)
%       y: predictant (3-D field, [lat,lon,time])
%       dof_reduction: in case the degrees of freedom is not (n-1, default) it can be specified 
%                      here as a scalar so [dof = n - dof_reduction]. Used for p-value.
%                      if empty --> dof_reduction = 1 by default
%
% output parameters (stored in a structure)
%       slope: regression slope (2-D field, [lat,lon])
%       intercept: regression intercept (2-D field, [lat,lon])
%       r2: explained variance (2-D field, [lat,lon])
%       p_val: p-value of regression (2-D field, [lat,lon])
%       rmse: rooted mean square error (2-D field, [lat,lon])
%       cor: correlation coefficient (2-D field, [lat,lon])
%       residuals: regression residuals (3-D field, [lat,lon,time])
%
% function [out] = Rean_regr(x,y,dof_reduction)
%
% FA 2015/11

function [out] = Rean_regr(x,y,dof_reduction)

%^^^^^^^^^^        check whether dof_reduction has been specified or use default      ^^^^^^^^^
if nargin<3
    dof_reduction = 1;
end

%^^^^^^^^^^^^^^^^^^^^^^           check if predictor has NaNs           ^^^^^^^^^^^^^^^^^^^^^^^
if sum(isnan(x))>0
    disp([num2str(sum(isnan(x))),' NaNs in the predictor will be ignored']);
end

%^^^^^^^^^^^^^^^^^^^^^^^          check if predictant has NaNs           ^^^^^^^^^^^^^^^^^^^^^^
if sum(sum(sum(isnan(y))))>0
    disp([num2str(sum(sum(sum(isnan(y))))),' NaNs in the field will be ignored']);
end

%^^^^^^^^^^^^^^^^^^^^^^^^^          internal                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
s_lat = size(y,1);
s_lon = size(y,2);
s_time = size(y,3);

%^^^^^^^^^^^^^^^^^^^^^            rearrage data into [time,lat x lon] space           ^^^^^^^^^
ys = shiftdim(y,2);
yr = reshape(ys,s_time,s_lat*s_lon);

%^^^^^^^^^^^^^^^^^^^^^       expand the predictor to matrix of same size             ^^^^^^^^^^
if size(x,1)<size(x,2)
    x = x';
end
if size(x,1)~=s_time
    error('variable and predictor must be of similar length');
end
xr = repmat(x,1,size(yr,2));

%^^^^^^^^^^^^^^^^^^^^^^           mean of each gridcell             ^^^^^^^^^^^^^^^^^^^^^^^^^^^
yrm = nanmean(yr);
xrm = nanmean(xr);

%^^^^^^^^^^^^^^^^^^^^^^            expand                           ^^^^^^^^^^^^^^^^^^^^^^^^^^^
xrme = repmat(xrm,s_time,1);
yrme = repmat(yrm,s_time,1);

%^^^^^^^^^^^^^^^^^^^^^             slope                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
slope = nansum((xr-xrme).*(yr-yrme))./nansum((xr-xrme).^2);

%^^^^^^^^^^^^^^^^^^^^^              intercept                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
intercept = yrm - slope.*xrm;

%^^^^^^^^^^^^^^^^^^^^^^             correlation              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
correlation = nansum((xr-xrme).*(yr-yrme))./sqrt(nansum((xr-xrme).^2).*nansum((yr-yrme).^2));

%^^^^^^^^^^^^^^^^^^^^^              r2                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
r2 = correlation.^2;

%^^^^^^^^^^^^^^^^^^^^    p-value based on  F-statistic (assuming n-1 dof)       ^^^^^^^^^^^^^^^
dof2 = s_time-dof_reduction;
dof1 = length(x)-dof2;
F = (r2./dof1)./((1-r2)./dof2);
p_value = 1-fcdf(F,dof1,dof2);

%^^^^^^^^^^^^^^^^^^^^               residuals & RMSE                ^^^^^^^^^^^^^^^^^^^^^^^^^^^
interc_exp = repmat(intercept,s_time,1);
slope_exp = repmat(slope,s_time,1);
residuals = yr - (interc_exp+xr.*slope_exp);
rmse = sqrt(nanmean(residuals.^2));

%^^^^^^^^^^^^^^^^^^^^^^            reshape to lat/lon               ^^^^^^^^^^^^^^^^^^^^^^^^^^^
slope_r = shiftdim(slope,1);
out.slope = reshape(slope_r,s_lat,s_lon);

intercept_r = shiftdim(intercept,1);
out.intercept = reshape(intercept_r,s_lat,s_lon);

rmse_r = shiftdim(rmse,1);
out.rmse = reshape(rmse_r,s_lat,s_lon);

p_value_r = shiftdim(p_value,1);
out.p_val = reshape(p_value_r,s_lat,s_lon);

correlation_r = shiftdim(correlation,1);
out.cor = reshape(correlation_r,s_lat,s_lon);

r2_r = shiftdim(r2,1);
out.r2 = reshape(r2_r,s_lat,s_lon);

residuals_r = shiftdim(residuals,1);
out.residuals = reshape(residuals_r,s_lat,s_lon,s_time);
end