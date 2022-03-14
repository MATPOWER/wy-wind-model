function model = wy_wind_build_model(wind_data)
%WY_WIND_BUILD_MODEL  Builds AR[1] model from wind speed data.
%
%   MODEL = WY_WIND_BUILD_MODEL(WIND_DATA)
%
%   Inputs:
%       WIND_DATA - (NW_ALL x NP_ALL) matrix of wind speeds (in m/s),
%           corresponding to NW_ALL specific sites, NP_ALL periods,
%           a particular NPD (number of periods per day), and a starting DT
%           specifying (year, month, day, period) (DT0)
%   Output:
%       MODEL - struct with fields:
%           ar1 - (NW_ALL x 1) vector of AR[1] coefficients for individual sites
%           ar1_total - scalar AR[1] coefficient for total wind

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Wooyoung Jeon
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

[nw_all, np_all] = size(wind_data);

%-----  WY your code goes here  -----
ar1_total = <scalar>
ar1 = <nw_all x 1 vector>




model = struct('ar1_total', ar1_total, 'ar1', ar1 );
