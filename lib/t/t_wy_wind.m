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

is_octave = exist('OCTAVE_VERSION', 'builtin') == 5;    %% Octave

if is_octave
    file_in_path_warn_id = 'Octave:data-file-in-path';
    s1 = warning('query', file_in_path_warn_id);
    warning('off', file_in_path_warn_id);
end

%% initialize parameters of interest
widx = [2;6;16];        %% wind sites of interest
nw = length(widx);      %% number of wind sites of interest
np = 12;                %% use a 12 period horizon
nb = 4;                 %% use 4 bins for forecasts
npd = 24;               %% number of periods per day in wind data
dt0 = [2004 1 1 1 0 0]; %% date of first period in wind data, 2004-01-01 1:00am
dt = [2004 8 1 0 0 0];  %% date of period of interest, 2004-08-01 12:00am
pidx1 = wy_wind_date2pidx(npd, dt0, dt);    %% scalar index of period of interest
init_rng(42);           %% initialize state of random number generator

t_begin(36+2*np, quiet);

%% load the historical data
t = 'wind_data_npcc';
s = load('wind_data_npcc'); %% 26303 x 16
wind_data = s.wind_data;
log_wind_data = log10(wind_data + 1);   %% convert raw wind to log(wind+1)
t_is(size(wind_data), [26303 16], 12, t);
t_ok(isfield(s, 'npd') && isequal(s.npd, 24), [t 'model.npd']);
t_ok(isfield(s, 'dt0') && isequal(s.dt0, [2004 1 1 1 0 0]), [t 'model.dt0']);

%% load the model
t = 'wind_model_npcc';
s = load('wind_model_npcc');
model = s.model;
t_ok(isstruct(model), [t 'model is struct']);
t_ok(isfield(model, 'type') && isequal(model.type, 1), [t 'model.type']);
t_ok(isfield(model, 'npd') && isequal(model.npd, 24), [t 'model.npd']);
t_ok(isfield(model, 'dt0') && isequal(model.dt0, [2004 1 1 1 0 0]), [t 'model.dt0']);
t_ok(isfield(model, 'ar1') && isequal(size(model.ar1), [16, 1]), [t 'model.ar1']);
t_ok(isfield(model, 'ols') && isequal(size(model.ols), [16, 9]), [t 'model.ols']);
t_ok(isfield(model, 'var_wnr') && isequal(size(model.var_wnr), [16, 16]), [t 'model.var_wnr']);
t_ok(isfield(model, 'ar1_total') && isequal(size(model.ar1_total), [1, 1]), [t 'model.ar1_total']);
t_ok(isfield(model, 'ols_total') && isequal(size(model.ols_total), [1, 9]), [t 'model.ols_total']);
t_ok(isfield(model, 'var_wnr_total') && isequal(size(model.var_wnr_total), [1, 1]), [t 'model.var_wnr_total']);

%% load the power curve
t = 'wind_power_curve_EIC.txt : ';
s2p = wy_wind_power_curve_data(5, 'wind_power_curve_EIC.txt');
t_is(size(s2p), [31, 2], 12, [t 'size']);
t_is(s2p(:, 1), [0:30]', 12, [t 's2p(:, 1)']);

%% create transition probabilities
tp = wy_wind_trans_probs(model, np, nb);

%% wind speed realization from data
wsr_data = wy_wind_realizations(log_wind_data, widx, pidx1, np);
wsr_data1 = wy_wind_realizations(log_wind_data, widx, pidx1+1, np);
wsr_data_c = (10.^wsr_data)-1;      %% convert log(wind+1) to wind speed

%% convert realization from speed to power
wpr_data = wy_wind_speed2power(wsr_data_c, 5);
wpr_data1 = wy_wind_speed2power(wsr_data1, s2p, 1);

%% wind speed realization from model
wsr_model = wy_wind_realizations(model, widx, pidx1, np);

%% convert realization from speed to power
wpr_model = wy_wind_speed2power(wsr_model, 5, 1);

%% generate forecast
ws0 = log_wind_data(pidx1, widx);
wsf = wy_wind_forecasts(model, widx, pidx1, ws0, np, nb);
wsf1 = wy_wind_forecasts(model, widx, pidx1+1, ws0, np, nb);
wsf1_c = (10.^wsf1)-1;  %% convert log(wind+1) to wind speed

%% convert forecast from speed to power
wpf = wy_wind_speed2power(wsf, s2p, 1);
wpf1 = wy_wind_speed2power(wsf1_c, 5);

%% load results
% if is_octave
%     oct = struct('wsr_model', wsr_model, 'wpr_model', wpr_model);
%     s = load('t_wy_wind_results');
%     ml = s.ml;
%     save -v7 t_wy_wind_results.mat tp wsr_data wpr_data ml oct wsf wpf wsr_data1 wpr_data1 wsf1 wpf1
% else                                        %% MATLAB
%     ml = struct('wsr_model', wsr_model, 'wpr_model', wpr_model);
%     save t_wy_wind_results.mat tp wsr_data wpr_data ml wsf wpf wsr_data1 wpr_data1 wsf1 wpf1
% end
s = load('t_wy_wind_results');

t = 'wy_wind_trans_probs : ';
t_is(length(tp), np, 12, [t 'length']);
for p = 1:np
    if p == 1
        t_is(size(tp{p}), [1 nb], 12, sprintf('size(tp{%d})', p));
    else
        t_is(size(tp{p}), [nb nb], 12, sprintf('size(tp{%d})', p));
    end
    t_is(tp{p}, s.tp{p}, 12, sprintf('%stp{%d}', t, p));
end

t = 'wy_wind_realizations(data, ...) : ';
t_is(size(wsr_data), [np, nw], 12, [t 'size(wsr)']);
t_is(wsr_data, s.wsr_data, 12, [t 'wsr_t0']);
t_is(wsr_data1, s.wsr_data1, 12, [t 'wsr_t1']);
t_is(wsr_data(2:end, :), wsr_data1(1:end-1, :), 12, [t 'wsr_t0(2:end,:) == wsr_t1(1:end-1,:)']);

t = 'wy_wind_speed2power : ';
t_is(size(wpr_data), [np, nw], 12, [t 'size(wpr)']);
t_is(wpr_data, s.wpr_data, 12, [t 'wpr_t0']);
t_is(wpr_data1, s.wpr_data1, 12, [t 'wpr_t1']);
t_is(wpr_data(2:end, :), wpr_data1(1:end-1, :), 12, [t 'wpr_t0(2:end,:) == wpr_t1(1:end-1,:)']);

t = 'wy_wind_realizations(model, ...) : ';
t_is(size(wsr_model), [np, nw], 12, [t 'size(wsr_model)']);
if is_octave
    t_is(wsr_model, s.oct.wsr_model, 12, [t 'wsr_model']);
else
    t_is(wsr_model, s.ml.wsr_model, 12, [t 'wsr_model']);
end

t = 'wy_wind_speed2power : ';
t_is(size(wsr_model), [np, nw], 12, [t 'size(wpr_model)']);
if is_octave
    t_is(wpr_model, s.oct.wpr_model, 12, [t 'wpr_model']);
else
    t_is(wpr_model, s.ml.wpr_model, 12, [t 'wpr_model']);
end

t = 'wy_wind_forecasts : ';
t_is(size(wsf), [np, nb, nw], 12, [t 'size(wsf)']);
t_is(wsf, s.wsf, 12, [t 'wsf_t0']);
t_is(wsf1, s.wsf1, 12, [t 'wsf_t1']);
t_ok(norm(wsf(2:end, :, 1) - wsf1(1:end-1, :, 1)) > 0.05, [t 'wsf_t0(2:end,:) ~= wsf_t1(1:end-1,:)']);

t = 'wy_wind_speed2power : ';
t_is(size(wpf), [np, nb, nw], 12, [t 'size(wpf)']);
t_is(wpf, s.wpf, 12, [t 'wpf_t0']);
t_is(wpf1, s.wpf1, 12, [t 'wpf_t1']);
t_ok(norm(wpf(2:end, :, 1) - wpf1(1:end-1, :, 1)) > 0.1, [t 'wpf_t0(2:end,:,1) ~= wpf_t1(1:end-1,:,1)']);

if is_octave
    warning(s1.state, file_in_path_warn_id);
end

t_end

function init_rng(N)
if exist('OCTAVE_VERSION', 'builtin') == 5  %% Octave
    randn('state', [1:N]);
else                                        %% MATLAB
    rng(N);
end
