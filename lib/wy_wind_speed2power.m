function wp = wy_wind_speed2power(ws, s2p, ws_type)
%WY_WIND_SPEED2POWER  Converts wind speeds to wind power
%
%   WP = WY_WIND_SPEED2POWER(WS, S2P)
%   WP = WY_WIND_SPEED2POWER(WS, IDX)
%   WP = WY_WIND_SPEED2POWER(WS, S2P, WS_TYPE)
%   WP = WY_WIND_SPEED2POWER(WS, IDX, WS_TYPE)
%
%   Inputs:
%       WS  - 2-D (NP x NW) or 3-D (NP x NB x NW) matrix of wind speeds
%       IDX - index of power curve to extract from default power curve data
%           Default IDX is 5, for multi-turbine
%       S2P - (M x 2) vector used as lookup table for converting wind
%           speed to wind power
%           col 1: wind speeds in m/s
%           col 2: corresponding power as fraction of installed capacity
%       WS_TYPE = how to interpret WS input
%           0 (default) - WS is provided as raw wind speed in m/s
%           1 - WS is provided as log(raw_wind_speed+1)
%
%   Output:
%       WP = 2-D or 3-D matrix of available wind power outputs as
%           fraction of installed capacity, same dimensions as WS

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Wooyoung Jeon & Ray Zimmerman
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

%% convert log(raw_wind_speed + 1) to raw_wind_speed, if indicated
if nargin == 3 && ws_type == 1
    ws = (10 .^ ws) - 1;
end

%% load power curve, if necessary
if isscalar(s2p)
    s2p = wy_wind_power_curve_data(s2p);
end

%% interpolate
wp = zeros(size(ws));
wp(:) = interp1(s2p(:, 1), s2p(:, 2), ws(:));
