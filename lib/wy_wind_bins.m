function [nb, bin_bounds] = wy_wind_bins(bins)
%WY_WIND_BINS  Returns NB and BIN_BOUNDS, given BINS
%
%   [NB, BIN_BOUNDS] = WY_WIND_BINS(BINS)
%
%   Input:
%       BINS  - bin specification, supplied as either:
%           (1) number of bins (NB), or
%           (2) (1 x NB-1) vector of bin boundaries (standard deviation
%               coefficients), where initial -Inf and final +Inf are assumed,
%               but not included
%
%   Output:
%       NB - number of bins
%       BIN_BOUNDS - 1 x (NB+1) bin boundaries, now including initial -Inf
%           and final +Inf

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Wooyoung Jeon
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

if isscalar(bins)
    nb = bins;

    %-----  WY your code goes here  -----
    % just need to assume something about the way the boundaries are
    % determined from a single number of bins NB
    bin_bounds = 
else
    nb = length(bins) + 1;
    bin_bounds = [-Inf bins Inf];
end
