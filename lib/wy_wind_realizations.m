function wsr = wy_wind_realizations(model, widx, pidx0, np)
%WY_WIND_REALIZATIONS  Returns wind speed realizations
%
%   WSR = WY_WIND_REALIZATIONS(MODEL, WIDX, PIDX0, NP)
%   WSR = WY_WIND_REALIZATIONS(WIND_DATA, WIDX, PIDX0, NP)
%
%   Inputs:
%       MODEL - struct with fields:
%           ar1 - (NW_ALL x 1) vector of AR[1] coefficients for individual sites
%           ar1_total - scalar AR[1] coefficient for total wind
%       WIND_DATA - (NW_ALL x NP_ALL) matrix of wind speeds (in m/s),
%           corresponding to NW_ALL specific sites, NP_ALL periods,
%           a particular NPD (number of periods per day), and a starting DT
%           specifying (year, month, day, period) (DT0)
%       WIDX  - (NW x 1) vector of indices of wind sites of interest
%       PIDX0 - scalar period index of first period of horizon of interest
%       NP    - number of periods of interest (e.g. for planning horizon)
%
%   Output:
%       WSR - (NW x NP) matrix of wind speed realizations

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Wooyoung Jeon
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

if ~isstruct(model)     %% extract realizations from data
    wsr = model(widx, pidx0:pidx0+np-1);
    return;
end

%% generate realizations from model
nw = length(widx);
wsr = zeros(nw, np);    %% initialize output matrix

%-----  WY your code goes here  -----
wsr(w, p) = 
