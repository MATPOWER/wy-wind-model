function tp = wy_wind_trans_probs(model, wind_data, ws0, np, bins)
%WY_WIND_TRANS_PROBS  Returns cell array of transition probabilities for MOST
%
%   TP = WY_WIND_TRANS_PROBS(MODEL, WS0, NP, BINS)
%
%   Inputs:
%       MODEL - struct with fields:
%           ar1 - (NW_ALL x 1) vector of AR[1] coefficients for individual sites
%           ar1_total - scalar AR[1] coefficient for total wind
%       WIND_DATA - (NW_ALL x NP_ALL) matrix of wind speeds (in m/s),
%           corresponding to NW_ALL specific sites, NP_ALL periods,
%           a particular NPD (number of periods per day), and a starting DT
%           specifying (year, month, day, period) (DT0)
%       WS0   - scalar, initial wind speed
%       NP    - number of periods of interest (e.g. for planning horizon)
%       BINS  - bin specification, supplied as either:
%           (1) number of bins (NB), or
%           (2) (1 x NB-1) vector of bin boundaries (standard deviation
%               coefficients), where initial -Inf and final +Inf are assumed,
%               but not included
%
%   Output:
%       TP - (1 x NP) cell array of transition probabilities, where
%           1st element is (1 x NB), rest are (NB x NB)

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Wooyoung Jeon
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

%% convert BIN spec argument to NB and BIN_BOUNDS
[nb, bin_bounds] = wy_wind_bins(bins);

%% initialize output array
tp = cell(1, np);

%-----  WY your code goes here  -----
tp{1} = 
...
tp{np} = 
