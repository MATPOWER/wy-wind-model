function wsf = wy_wind_forecasts(model, widx, pidx0, ws0, np, bins)
%WY_WIND_FORECASTS  Returns wind speed forecast bin means
%
%   WSF = WY_WIND_FORECASTS(MODEL, WIDX, PIDX0, WS0, NP, BINS)
%
%   Inputs:
%       MODEL - struct with fields:
%           ar1 - (NW_ALL x 1) vector of AR[1] coefficients for individual sites
%           ols - (NW_ALL x 9) matrix of OLS estimation parameters for
%               individual wind sites: [C CY1 SY1 CY2 SY2 CD1 SD1 CD2 SD2]
%           var_wnr - (NW_ALL x NW_ALL) covariance matrix for individual sites
%           ar1_total - scalar AR[1] coefficient for total wind
%           ols_total - 1 x 9 vector of OLS estimation parameters for total wind
%           var_wnr - scalar variance for total wind
%       WIDX  - (NW x 1) vector of indices of wind sites of interest
%       WS0 - (NW x 1) vector of initial wind, units must be consistent with
%           those used by MODEL
%       PIDX0 - scalar period index of first period of horizon of interest
%       NP    - number of periods of interest (e.g. for planning horizon)
%       BINS  - bin specification, supplied as either:
%           (1) number of bins (NB), or
%           (2) (1 x NB-1) vector of bin boundaries (standard deviation
%               coefficients), where initial -Inf and final +Inf are assumed,
%               but not included
%   Output:
%       WSF - (NP x NB x NW) 3D matrix of bin means of wind forecast, units are
%           consistent with those used by MODEL

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Wooyoung Jeon
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

%% convert BIN spec argument to NB and BIN_BOUNDS
[nb, bin_bounds] = wy_wind_bins(bins);

%% generate forecasts from model
nw = length(widx);
wsf = zeros(np, nb, nw);    %% initialize output matrix

probcum = normcdf_nst(bin_bounds);  % cdf for each sd
probref = probcum(1:end-1) + diff(normcdf_nst(bin_bounds))/2;   % find mean prob location for each bin in normal dist
bin_mean = norminv_nst(probref);    % find mean location for each bin in normal dist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. initial setup for dataset and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply estimates for total wind
coef_cycle = model.ols(:,2:end);    % coefficient for cycles, (8x1)
coef_mean = model.ols(:,1);         % coefficient for mean(constant), (1x1)

rho = model.ar1;
sd_wnr= sqrt(model.var_wnr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. inputs from econometric model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cycle information defined
% PERIODS
PY1 = 8766;
PY2 = PY1 / 2;
PD1 = 24;
PD2 = PD1 / 2;

% adjustment for calender cycle input, push 1 hour to the next
% IMPORTANT: cycle hour, 1hour shifted as estimation is done this way
% hour 0 to hour 24, 25hours,
% hour 0 needed for ar(1) process. at(t-1) is needed
shift = 1;
tt2=[pidx0+shift:1:pidx0+np+shift]';

% cosine and sine of full year, half year, full day, half day
c_y1 = cos( (2*pi()/ PY1) * tt2 );
s_y1 = sin( (2*pi()/ PY1) * tt2 );
c_y2 = cos( (2*pi()/ PY2) * tt2 );
s_y2 = sin( (2*pi()/ PY2) * tt2 );
c_d1 = cos( (2*pi()/ PD1) * tt2 );
s_d1 = sin( (2*pi()/ PD1) * tt2 );
c_d2 = cos( (2*pi()/ PD2) * tt2 );
s_d2 = sin( (2*pi()/ PD2) * tt2 );

% cycle variables in matrix
var_cycle = [c_y1, s_y1, c_y2, s_y2, c_d1, s_d1, c_d2, s_d2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. generates forecasted bin, ff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% yy:  yact : realized LWIND
% ww: yfit : ols fitted LWIND => MEAN
% uu: yres : ols residual of LWIND
% yres_fit : ar(1) fitted LWIND
% ee: wnr : ar(1) white noise residual of LWIND
% ff: yfor : forecasted LWIND = yfit+yres_fit = ww + yres_fit

for i=1:nw
    yINIT = ws0(i);
    ww= coef_mean(widx(i)) + var_cycle * coef_cycle(widx(i),:)';    %yfit
    uINIT = yINIT - ww(1);

    sumrho = (rho(widx(i)).^2).^[0:np-1]';

    for t=1:np
        var_for(t,1) = sd_wnr(widx(i), widx(i))^2 * sum(sumrho(1:t,1));
    end

    sd_for = sqrt(var_for); % sqrt(var[forecasted LWIND])

    % forecasted LWIND
    ff0 = zeros(np+1,1);
    ff0(1) = yINIT;
    for t=1:np
        ff0(t+1) = ww(t+1) - rho(widx(i))*ww(t) + rho(widx(i))*ff0(t);
    end
    ff = ff0(2:end);

    bm = ff0(2:end) + sd_for * bin_mean;

    % uuu : AR residual for bin mean i in hour t-1
    % bf : forecast for bin mean i in hour t

    uuu0 = zeros(np+1,nb);  % t starting at 0, 25 hours

    midx = floor((nb+1)/2); % index for central bin

    uuu0(1,midx) = uINIT;
    for t=2:np+1
        uuu0(t,:) = bm(t-1,:) - ff0(t,1);
    end

    uuu = uuu0(2:end,:);

    bf0 = zeros(np-1, nb);
    for t=2:np
        bf0(t-1,:) = ff(t,1) + rho(widx(i))*uuu0(t,:);
    end

    bf = [ones(1,nb) * ff(1); bf0];

    bf_all(:,:,i) = bf;
    wsf(:,:,i) = bf;
end
