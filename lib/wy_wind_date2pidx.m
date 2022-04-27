function pidx = wy_wind_date2pidx(npd, dt0, dt)
%WY_WIND_DATE2PIDX  Converts a date vector to a raw period index
%
%   PIDX = WY_WIND_DATE2PIDX(NPD, DT0, DT)
%
%   Inputs:
%       NPD - number of periods per day (typically 24, for hourly data)
%       DT0 - 1 x 6 date vector specifying initial period of raw historical data
%           (same format as DT)
%       DT = 1 x 6 date vector (YR, MO, DAY, HR, MIN, SEC), specifying a
%           specific period in the raw historical wind data, where
%           YR - 4-digit year
%           MO - month-of-year (1-12)
%           DAY - day-of-month (1-31)
%           MIN - minutes (0-59)
%           SEC - seconds (0-59)
%
%   Output:
%       PIDX = scalar period index for raw historical wind data

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Ray Zimmerman
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

pidx = round(npd * (datenum(dt) - datenum(dt0))) + 1;
