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
    switch nb
        case 1
            bin_bounds = [-inf inf];
        case 2
            bin_bounds = [-inf 0 inf];
        case 3
            bin_bounds = [-inf -1 1 inf];
        case 4
            bin_bounds = [-inf -1 0 1 inf];
        case 5
            bin_bounds = [-inf -2 -1 1 2 inf];
        case 6
            bin_bounds = [-inf -2 -1 0 1 2 inf];
        case 7
            bin_bounds = [-inf -3 -2 -1  1 2 3 inf];
        case 8
            bin_bounds = [-inf -3 -2 -1 0 1 2 3 inf];
        case 9
            bin_bounds = [-inf -4 -3 -2 -1 1 2 3 4 inf];
        case 10
            bin_bounds = [-inf -4 -3 -2 -1 0 1 2 3 4 inf];
        otherwise
            error('wy_wind_bins: BINS = %d, scalar BINS must be between 1 and 10', nb);
    end
else
    nb = length(bins) + 1;
    bin_bounds = [-Inf bins Inf];
end
