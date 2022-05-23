function wsr = wy_wind_realizations(model, widx, pidx0, np)
%WY_WIND_REALIZATIONS  Returns wind speed realizations
%
%   WSR = WY_WIND_REALIZATIONS(MODEL, WIDX, PIDX0, NP)
%   WSR = WY_WIND_REALIZATIONS(WIND_DATA, WIDX, PIDX0, NP)
%
%   If the first argument is a struct (MODEL), a new realization is genereted
%   using the model. On the other hand, if the first argument is a matrix
%   (WIND_DATA) the realization is extracted from the historical data.
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
%       WIND_DATA - (NP_ALL x NW_ALL) matrix of wind speeds, corresponding
%           to NP_ALL periods, NW_ALL specific sites, a particular NPD (number
%           of periods per day), and a starting date/time (DT0)
%       WIDX  - (NW x 1) vector of indices of wind sites of interest
%       PIDX0 - scalar period index of first period of horizon of interest
%       NP    - number of periods of interest (e.g. for planning horizon)
%
%   Output:
%       WSR - (NP x NW) matrix of wind speed realizations, the units are
%           consistent with those used by MODEL or WIND_DATA, respectively.

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Wooyoung Jeon
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

if isnumeric(model) %% extract realization from data
    wsr = model(pidx0:pidx0+np-1, widx);
else                %% generate realization from model
    nw = length(widx);

    % recreate var_wnr for selected wind site, widx
    temp_var = model.var_wnr(widx, widx);

    % copy lower triangular part to upper triangular part to make it symmatric
    % var-covar matrix
    var_wnr=tril(temp_var,-1)'+temp_var;

    % generate randomized wnr based on normal distribution using var-covar
    % matrix
    gen_wnr = mvnrnd_nst(zeros(1,nw),var_wnr,np);

    % ar(1) part
    % for t=1:np
    %   ar_wnr(t,:) = gen_wnr(t,:) .* model.ar1^t;
    % end

    temp1 = zeros(np, np, nw);
    ar_sum = zeros(np, nw);
    for t=1:nw
        for i=1:np
            for j=1:np
                % create matrix of [e1 e1*ar1 e1*ar1^2 ... e1*ar1^23] for each
                % e_i
                temp1(i,j,t) = gen_wnr(i,t) .* model.ar1(widx(t))^(j-1);
            end
        end
        temp2=squeeze(temp1(:,:,t))';

        % flip the matrix and sum diagonal of each of increasing matrix size
        for i=1:np
            temp3=flip(squeeze(temp2(1:i,1:i)));
            ar_sum(i,t)=sum(diag(temp3));   % (np x nw)
        end
    end

    % ols part : computing mean logwind
    % cycle information defined
    % PERIODS
    PY1 = 8766;
    PY2 = PY1 / 2;
    PD1 = 24;
    PD2 = PD1 / 2;

    % adjustment for calender cycle input, push 1 hour to the next
    % IMPORTANT: cycle hour, 1hour shifted as estimation is done this way
    % hour 1 to hour 24, 24hours,
    % hour 0 needed for ar(1) process. at(t-1) is needed
    shift = 1;
    tt2=[pidx0+shift+1:1:pidx0+np+shift]';

    % cosine and sine of full year, half year, full day, half day
    c_y1 = cos( (2*pi()/ PY1) * tt2 );
    s_y1 = sin( (2*pi()/ PY1) * tt2 );
    c_y2 = cos( (2*pi()/ PY2) * tt2 );
    s_y2 = sin( (2*pi()/ PY2) * tt2 );
    c_d1 = cos( (2*pi()/ PD1) * tt2 );
    s_d1 = sin( (2*pi()/ PD1) * tt2 );
    c_d2 = cos( (2*pi()/ PD2) * tt2 );
    s_d2 = sin( (2*pi()/ PD2) * tt2 );

    % cycle variables in matrix
    var_cycle = [c_y1, s_y1, c_y2, s_y2, c_d1, s_d1, c_d2, s_d2];

    mean_wind = zeros(np, nw);
    for i=1:nw
        mean_wind(:,i)= model.ols(widx(i),1) + var_cycle * model.ols(widx(i),2:end)'; %yfit, (np x nw)
    end

    % ols + ar(1)
    wsr = mean_wind + ar_sum;
end
