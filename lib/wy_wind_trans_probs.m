function tp = wy_wind_trans_probs(model, np, bins)
%WY_WIND_TRANS_PROBS  Returns cell array of transition probabilities for MOST
%
%   TP = WY_WIND_TRANS_PROBS(MODEL, WS0, NP, BINS)
%
%   Generates transition probabilities for MOST from a time series model.
%
%   Inputs:
%       MODEL - struct with fields:
%           type -  type of wind speed model
%               0 = based on raw_wind_speed in m/s
%               1 = based on log10(raw_wind_speed + 1)
%           npd - number of periods per day (24 for hourly model)
%           dt0 - Matlab date vector corresponding to first period of data
%               from which model was created
%           ar1 - (NW_ALL x 1) vector of AR[1] coefficients for individual sites
%           ols - (NW_ALL x 9) matrix of OLS estimation parameters for
%               individual wind sites: [C CY1 SY1 CY2 SY2 CD1 SD1 CD2 SD2]
%           var_wnr - (NW_ALL x NW_ALL) covariance matrix for individual sites
%           ar1_total - scalar AR[1] coefficient for total wind
%           ols_total - 1 x 9 vector of OLS estimation parameters for total wind
%           var_wnr - scalar variance for total wind
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
[nb, sb] = wy_wind_bins(bins);

rho=model.ar1_total;                % ar1 coefficient of total model
sd_wnr = sqrt(model.var_wnr_total); % sd of wnr of total model

sumrho = (rho^2).^[0:np-1]';
vf = zeros(np, 1);
for t=1:np
    vf(t,1) = sd_wnr^2 * sum(sumrho(1:t,1));
end

sdf = sqrt(vf);
sm = zeros(nb, 1);
for i = 1:nb
    sm(i)=norminv_nst(normcdf_nst(sb(i)) + (normcdf_nst(sb(i+1)) - normcdf_nst(sb(i))) / 2);
end

midx = floor((nb+1)/2);     % index for central bin

% tp(1) : probability of being in each bin
tp1 = zeros(nb, nb, np);
for j=1:nb
    tp1(midx,j,1)= normcdf_nst(sb(j+1)) - normcdf_nst(sb(j));
end
for t=2:np
    for i=1:nb
        for j=1:nb
            tp1(i,j,t) = normcdf_nst(sb(j+1)*sqrt(vf(t))-rho*sm(i)*sqrt(vf(t-1)),0,sdf(1,1))- normcdf_nst(sb(j)*sqrt(vf(t))-rho*sm(i)*sqrt(vf(t-1)),0,sdf(1,1));
        end
    end
end

% tp : cell{1 x np}
% cell(1) : {1 x nb}, cell(2:np) : {nb x nb}
tp = cell(1,np);
tp{1} = squeeze(tp1(midx,:,1))';
for t=2:np
    tp{t} = squeeze(tp1(:,:,t))';
end
