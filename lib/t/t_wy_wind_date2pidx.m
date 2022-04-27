function t_wy_wind_date2pidx(quiet)
%T_WY_WIND_DATE2PIDX  Tests WY-Wind-Model WY_WIND_DATE2PIDX function.

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Ray Zimmerman
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

if nargin < 1
    quiet = 0;
end

t_begin(8, quiet);

t_is(wy_wind_date2pidx(24, [2000 12 25 12], [2000 12 25 12]), 1, 12, '1st index');
t_is(wy_wind_date2pidx(24, [2020 2 28 0], [2020 2 28 0]), 1, 12, '1st index');
t_is(wy_wind_date2pidx(24, [2020 2 28 4], [2020 3 1 4]), 49, 12, 'leap year');
t_is(wy_wind_date2pidx(24, [2020 2 28 4], [2021 3 1 4]), 25+366*24, 12, 'leap year');
t_is(wy_wind_date2pidx(24, [2021 2 28 6], [2021 3 1 6]), 25, 12, 'non-leap year');
t_is(wy_wind_date2pidx(12, [2021 2 28 0], [2021 3 1 0]), 13, 12, 'non-leap year');
t_is(wy_wind_date2pidx(3, [2007 6 15 2], [2021 7 16 0]), (365*14+30+4)*3+2, 12, 'multiple years');
t_is(wy_wind_date2pidx(24, [2004 1 1 1], [2004 8 1 0]), 5112, 12, 'WY example');

t_end;
