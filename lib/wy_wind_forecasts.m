function wsf = wy_wind_forecasts(model, widx, pidx0, np, bins)
%WY_WIND_FORECASTS  Returns wind speed forecast bin means
%
%   WSF = WY_WIND_FORECASTS(MODEL, WIDX, PIDX0, NP, BINS)
%
%   Inputs:
%       MODEL - struct with fields:
%           ar1 - (NW_ALL x 1) vector of AR[1] coefficients for individual sites
%           ar1_total - scalar AR[1] coefficient for total wind
%       WIDX  - (NW x 1) vector of indices of wind sites of interest
%       PIDX0 - scalar period index of first period of horizon of interest
%       NP    - number of periods of interest (e.g. for planning horizon)
%       BINS  - bin specification, supplied as either:
%           (1) number of bins (NB), or
%           (2) (1 x NB-1) vector of bin boundaries (standard deviation
%               coefficients), where initial -Inf and final +Inf are assumed,
%               but not included
%   Output:
%       WSF - (NW x NP x NB) 3D matrix of bin means of wind speed forecast

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
wsf = zeros(nw, np, nb);    %% initialize output matrix

%-----  WY your code goes here  -----
wsf(w, p, b) = 
