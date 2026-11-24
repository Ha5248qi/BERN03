%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%             PCA analysis on Reanalysis data             %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script performs PCA analysis on gridded reanalysis data
% The data will be weighted by latitude prior
% only gridcells that do not contain a single NaN in their timeseries will be considered
%
% Input Parameters
%   lat:        Latitude Vector
%   lon:        Longitude Vector
%   data:       3-D Data Matrix, shaped (lat,lon,time)
%   no_of_pcs:  number of PCs to be calculated (default=10)
%   option:     'detrend' (for linear detrending before PCA)
%
% Output parameters
%   COEFF:      "loadings" for each gridcell and PC, shaped (lat,lon,PC)
%   PC:         time series of PCs, shaped (time,PC)
%   VARexp:     explained variance for each PC
%   LAMBDA:     Eigenvectors
%   d_LAMBDA:   Sampling error according to North (1982)
%
% FA 2015/07/28 (largely based on Jesper Sjoltes eof_global_TCR.mat)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [COEFF,PC,VARexp,LAMBDA,d_LAMBDA] = Rean_PCA_anomalies(lat,lon,data,no_of_pcs,option)

% check input arguments
if nargin==3
    no_of_pcs = 10;
    option = 'none';
end
if nargin == 4
    option = 'none';
end
if nargin == 5 && ~strcmpi(option,'detrend')
    error('Unknown Option.');
end

%^^^^^^^^^^^^^^^^^            get sizes                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
x_size = size(data,1);
y_size = size(data,2);
t_size = size(data,3);

%^^^^^           check if sizes are consistent with data = [lat,lon,time]    ^^^^^^^^^^^^^^^^^^
if x_size~=length(lat) || y_size~=length(lon) || t_size ~= t_size
    error('Input data must be in [lat,lon,time] space.');
end

%^^^^^^^^^^^^^^^^^            make lat-lon grid           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
lat_mat = double(meshgrid(lat,lon)');
lon_mat = double(meshgrid(lon,lat));

%^^^^^^^^^^^^^^^^^^           expand lat_mat to time space           ^^^^^^^^^^^^^^^^^^^^^^^^^^
lat_mat_t = repmat(lat_mat,[1,1,t_size]);
lon_mat_t = repmat(lon_mat,[1,1,t_size]);
%^^^^^^^^^^^^^^^^          NaN index                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
nan_mat = isnan(data);

%^^^^^^^^^^^^^^^^^^^              reshape data to [time,lat*lon]        ^^^^^^^^^^^^^^^^^^^^^^^
data_s = shiftdim(data,2);
data_r = reshape(data_s,[t_size,x_size*y_size]);

lat_mat_s = shiftdim(lat_mat_t,2);
lat_mat_r = reshape(lat_mat_s,[t_size,x_size*y_size]);

lon_mat_s = shiftdim(lon_mat_t,2);
lon_mat_r = reshape(lon_mat_s,[t_size,x_size*y_size]);

nan_mat_s = shiftdim(nan_mat,2);
nan_mat_r = reshape(nan_mat_s,[t_size,x_size*y_size]);

%^^^^^^^^^^^^^^^^^^^^           get grid cells that do not contain nans      ^^^^^^^^^^^^^^^^^^
nan_cell = sum(nan_mat_r)==0;

%^^^^^^^^^^^^^^^^^^^^            detrend                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
if strcmpi(option,'detrend')
    data_r = detrend(data_r);
end

%^^^^^^^^^^^^^^^^^^^             zscores              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%data_r = zscore(data_r);
data_r = bsxfun(@minus, data_r,mean(data_r));

%^^^^^^^^^^^^^^^^^^^^          latitude weighting          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
data_r = data_r.*sqrt(cos(lat_mat_r*pi/180));

%^^^^^^^^^^^^^^^^^^^          take cells that do not contain NaNs       ^^^^^^^^^^^^^^^^^^^^^^^
data_rn = data_r(:,nan_cell);

%^^^^^^^^^^^^^^^^^^^^          covariance matrix             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Rcov = data_rn'*data_rn;

%^^^^^^^^^^^^^^               Eigenvalues                                        ^^^^^^^^^^^^^^
Rcov = double(Rcov); 

k = no_of_pcs; % number of eigenvectors to calculate
[COEFF,LAMBDA] = eigs(Rcov,k);

%^^^^^^^^^^^^^^^                   PC time series                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Find expansion coefficients corresponding to eigenvalues number i
PC = data_rn*COEFF;

%^^^^^^^^^^^^^^^^               Explained variance of PCs            ^^^^^^^^^^^^^^^^^^^^^^^^^^
VARexp = diag(LAMBDA)/trace(LAMBDA);

%^^^^^^^^^^^^^^^^^^               User output                 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
COEFF = COEFF';
COEFF_out = NaN(no_of_pcs,size(data_r,2));
COEFF_out(:,nan_cell) = COEFF;

COEFF_out = shiftdim(COEFF_out,1);
COEFF_out = reshape(COEFF_out,[x_size, y_size, no_of_pcs]);
COEFF = COEFF_out;

N = t_size; % Numberof years in sample
E = VARexp; % Relative size of eigenvector
d_LAMBDA = E*sqrt(2/N); % 'Sampling error' according to North (1982)
end