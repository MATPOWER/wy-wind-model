function s2p = wy_wind_power_curve_data(idx, fname)
%WY_WIND_POWER_CURVE_DATA  Converts wind speeds to wind power
%
%   S2P = WY_WIND_POWER_CURVE_DATA(IDX)
%   S2P = WY_WIND_POWER_CURVE_DATA(IDX, FNAME)
%
%   Inputs:
%       IDX - index of power curve to extract from power curve data
%           Default IDX is 5, for multi-turbine
%       FNAME - (optional) name of raw text file containing power curve data
%           as a matrix where the first column is wind speed in m/s, and
%           subsequent columns (2..N+1) contain the power curves corresponding
%           to IDX (1..N) as fraction of installed capacity
%           Default file name is 'WindPowerCurveIEC.txt'
%               col 1 is wind speeds from 0 to 30 m/s
%               cols 2-6 contain 5 power curves
%                   1: IEC1, 2:IEC2, 3:IEC3, 4:Offshore, 5:Multi-turbine
%
%   Output:
%       S2P - ((WS_MAX+1) x 2) vector used as lookup table for converting wind
%           speed to wind power
%           col 1: wind speeds in m/s
%           col 2: corresponding power as fraction of installed capacity

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Ray Zimmerman
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

%% default inputs
if nargin < 2
    fname = 'WindPowerCurveIEC.txt';
    if nargin < 1
        idx = 5;
    end
end

data = load(fname);
s2p = [data(:, 1) data(:, idx+1)];
