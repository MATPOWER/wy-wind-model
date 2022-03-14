function wp = wy_wind_speed2power(ws, s2p)
%WY_WIND_SPEED2POWER  Converts wind speeds to wind power
%
%   WP = WY_WIND_SPEED2POWER(WS, S2P)
%
%   Inputs:
%       WS  - (NW x NP x NB) 3D matrix of wind speeds
%       S2P - (1 x (WS_MAX+1)) vector (to be used for all sites of interest) or
%           (NW x (WS_MAX+1)) matrix of fractions (for NW individual sites),
%           used as lookup table for converting wind speed to wind power, as
%           fraction of installed capacity, where 1st and last cols correspond
%           to wind speeds of 0 m/s and ws_max m/s, respectively
%   Output:
%       WP = (NW x NP x NB) 3D matrix of available wind power outputs as
%           fraction of installed capacity

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Wooyoung Jeon
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

%-----  WY your code goes here  -----
