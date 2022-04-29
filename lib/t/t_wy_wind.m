function t_wy_wind(quiet)
%T_WY_WIND  Tests WY-Wind-Model functions.

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

if have_feature('octave')
    file_in_path_warn_id = 'Octave:data-file-in-path';
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end

%% initialize parameters of interest
widx = [2;6;16];        %% wind sites of interest
np = 12;                %% use a 12 period horizon
bins = 4;               %% use 4 bins for forecasts
npd = 24;               %% number of periods per day in wind data
dt0 = [2004 1 1 1 0 0]; %% date of first period in wind data, 2004-01-01 1:00am
dt = [2004 8 1 0 0 0];  %% date of period of interest, 2004-08-01 12:00am
pidx0 = wy_wind_date2pidx(npd, dt0, dt);    %% scalar index of period of interest

t_begin(3+np, quiet);

%% load the historical data
% wind_data = load('winddata_npcc');  %% 26303 x 16
s = load('winddata_npcc');  %% 26303 x 16
wind_data = s.wind_data;

%% load the model
% model = load('model_npcc');
s = load('model_npcc');
model = s.model;

%% create transition probabilities
tp = wy_wind_trans_probs(model, np, bins);

%% wind speed realization from data
wsr_data = wy_wind_realizations(wind_data, widx, pidx0, np);

%% convert realization from speed to power
s2p = load('WindPowerCurveIEC.txt');
s2p = s2p(:, 6)';       %% select curve for multi-turbine
% wpr = wy_wind_speed2power(wsr_model, s2p);
wpr_data = wy_wind_speed2power(wsr_data, 5);

% %% wind speed realization from model
% wsr_model = wy_wind_realizations(model, widx, pidx0, np);
% 
% %% convert realization from speed to power
% % wpr = wy_wind_speed2power(wsr_model, s2p);
% wpr = wy_wind_speed2power(wsr_model, 5);

% %% generate forecast
% % wsf = wy_wind_forecasts(model, widx, pidx0, np, bins);
% wsf = wy_wind_forecasts(wind_data, model, widx, pidx0, np, bins);
% 
% %% convert forecast from speed to power
% % wpf = wy_wind_speed2power(wsf, s2p);
% wpf = wy_wind_speed2power(wsf, 5);

%% load results
s = load('t_wy_wind_results');

t = 'wy_wind_trans_probs : ';
t_is(length(tp), np, 12, [t 'length']);
for p = 1:np
    t_is(tp{p}, s.tp{p}, 12, sprintf('%stp{%d}', t, p));
end

t = 'wy_wind_realizations(data, ...)';
t_is(wsr_data, s.wsr_data, 12, t);

t = 'wy_wind_speed2power';
t_is(wpr_data, s.wpr_data, 12, t);


if have_feature('octave')
    warning(s1.state, file_in_path_warn_id);
end

t_end;
