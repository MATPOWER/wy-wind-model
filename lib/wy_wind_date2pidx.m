function pidx = wy_wind_date2pidx(npd, dt0, dt)
%WY_WIND_DATE2PIDX  Converts a date vector to a raw period index
%
%   PIDX = WY_WIND_DATE2PIDX(NPD, DT0, DT)
%
%   Inputs:
%       NPD - number of periods per day (typically 24, for hourly data)
%       DT0 - 1 x 4 vector specifying initial period of raw historical data
%           (same format as DT)
%       DT = 1 x 4 vector (YR, MO, DAY, P), specifying a specific period in
%           the raw historical wind data, where
%           YR - 4-digit year
%           MO - month-of-year (1-12)
%           DAY - day-of-month (1-31)
%           P - period of day (0-(NPD-1)) (e.g. 0-23 for hourly, i.e. NPD=24)
%
%   Output:
%       PIDX = scalar period index for raw historical wind data

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Wooyoung Jeon
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

%-----  WY your code goes here  -----
pidx = 