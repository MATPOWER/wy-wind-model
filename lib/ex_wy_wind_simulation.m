%EX_WY_WIND_SIMULATION Example of wind inputs for receding horizon

%% load historical wind data
s = load('wind_data_npcc');
wind_data = s.wind_data;
log_wind_data = log10(wind_data + 1);

widx = [1;3;5;7];   %% wind sites of interest
np = 24;            %% num of periods
nb = 5;             %% num of bins
pidx00 = wy_wind_date2pidx(24, [2004 1 1 1 0 0], [2004 7 27 8 0 0]);    %% 5000
nw = length(widx);  %% num of wind sites of interest
nrh = 24;           %% length of receding horizon window

%% create wind model
wm = wy_wind_model('wind_model_npcc', widx);

%% initialize outputs
tp = cell(nrh, 1);
wsr = cell(nrh, 1);
wpr = cell(nrh, 1);
wsf = cell(nrh, 1);
wpf = cell(nrh, 1);

%% generate wind inputs for defined receding horizon window
for t = 1:nrh
    %% create transition probability matrices
    tp{t} = wm.transition_probs(np, nb);

    %% starting period index for window t
    pidx0 = pidx00 - 1 + t;

    %% generate realized wind power
    wsr{t} = wm.realizations(pidx0, np);
    wpr{t} = wm.speed2power(wsr{t});

    %% generate forecasted wind power
    ws0 = log_wind_data(pidx0-1, widx);
    wsf{t} = wm.forecasts(pidx0, ws0, np, nb);
    wpf{t} = wm.speed2power(wsf{t});
end

%% plot log forecasted wind for first window
for i = 1:nw
    figure(i)
    plot(wsf{1}(:, :, i));
end
