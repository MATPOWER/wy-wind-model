classdef wy_wind_model < handle
%WY_WIND_MODEL  WY-Wind-Model class
%
%   WM = WY_WIND_MODEL(MODEL_FNAME)
%   WM = WY_WIND_MODEL(MODEL_FNAME, WIDX)
%   WM = WY_WIND_MODEL(MODEL_FNAME, WIDX, S2P)
%
%   Properties
%       type          - type of wind speed model
%                          0 - based on raw wind speed in m/s
%                          1 - based on log10(raw_wind_speed + 1)
%       nw            - number of wind sites
%       widx          - indices of wind sites in original data (use ':' for all)
%       ar1           - nw x 1 vector of AR(1) coefficients for indiv sites
%       ols           - nw x 9 matrix of OLS estimation parameters for
%                       individual sites: [C CY1 SY1 CY2 SY2 CD1 SD1 CD2 SD2]
%       var_wnr       - (nw x nw) covariance matrix for individual sites
%       ar1_total     - scalar AR[1] coefficient for total wind
%       ols_total     - 1 x 9 vector of OLS estimation parameters for total wind
%       var_wnr_total - variance for total wind
%       s2p           - power curve table for converting wind speed to power
%
%   Methods
%       wy_wind_model()
%           wm = wy_wind_model(model_fname)
%           wm = wy_wind_model(model_fname, widx)
%           wm = wy_wind_model(model_fname, widx, s2p)
%       transition_probs()
%           tp = wm.transition_probs(np, bins)
%       realizations()
%           wsr = wm.realizations(pidx0, np)
%           wsr = wm.realizations(pidx0, np, wind_data)
%       forecasts()
%           wsf = wm.forecasts(pidx0, ws0, np, bins)
%       speed2power()
%           wp = wm.speed2power(ws)
%       display() - called to display object on command line

%   WY-Wind-Model
%   Copyright (c) 2022, Wooyoung Jeon, Ray Zimmerman
%   by Wooyoung Jeon & Ray Zimmerman
%
%   This file is part of WY-Wind-Model.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/wy-wind-model for more info.

    properties
        type = 0;       %% type of wind speed model
                        %%    0 - based on raw wind speed in m/s
                        %%    1 - based on log10(raw_wind_speed + 1)
        npd             %% number of periods per day (24 for hourly model)
        dt0             %% Matlab date vector corresponding to first period of
                        %% data from which model was created
        nw              %% number of wind sites
        widx            %% indices of wind sites in original data
        ar1             %% nw x 1 vector of AR(1) coefficients for indiv sites
        ols             %% nw x 9 matrix of OLS estimation parameters for
                        %% individual sites: [C CY1 SY1 CY2 SY2 CD1 SD1 CD2 SD2]
        var_wnr         %% (nw x nw) covariance matrix for individual sites
        ar1_total       %% scalar AR[1] coefficient for total wind
        ols_total       %% 1 x 9 vector of OLS estimation parameters for total wind
        var_wnr_total   %% variance for total wind
        s2p             %% power curve table for converting wind speed to power
    end

    methods
        %% constructor
        function obj = wy_wind_model(model_fname, widx, s2p)
            %% load model from file
            s = load(model_fname);
            model = s.model;

            %% select wind sites of interest
            if nargin < 2 || (ischar(widx) && strcmp(widx, ':'))
                obj.nw = length(model.ar1);
                widx = [1:obj.nw]';
            else
                obj.nw = length(widx);
            end
            obj.widx = widx;
            model.ar1 = model.ar1(widx, :);
            model.ols = model.ols(widx, :);
            model.var_wnr = model.var_wnr(widx, widx);

            %% assign model parameters
            obj.type = model.type;
            obj.npd = model.npd;
            obj.dt0 = model.dt0;
            obj.ar1 = model.ar1;
            obj.ols = model.ols;
            obj.var_wnr = model.var_wnr;
            obj.ar1_total = model.ar1_total;
            obj.ols_total = model.ols_total;
            obj.var_wnr_total = model.var_wnr_total;
        end

        function display(obj)
            fprintf('WY-Wind-Model for %d sites using ', obj.nw);
            switch obj.type
                case 0
                    fprintf('raw wind speed (m/s)\n');
                case 1
                    fprintf('log10(raw_wind_speed (m/s) + 1)\n');
                otherwise
                    fprintf('<UNKNOWN TYPE>\n');
            end
        end

        function tp = transition_probs(obj, np, bins)
            tp = wy_wind_trans_probs(obj, np, bins);
        end

        function wsr = realizations(obj, pidx0, np, wind_data)
            if nargin > 3 && ~isempty(wind_data)
                wsr = wy_wind_realizations(wind_data, obj.widx, pidx0, np);
            else
                wsr = wy_wind_realizations(obj, [1:obj.nw]', pidx0, np);
            end
        end

        function wsf = forecasts(obj, pidx0, ws0, np, bins)
            wsf = wy_wind_forecasts(obj, [1:obj.nw]', pidx0, ws0, np, bins);
        end

        function wp = speed2power(obj, ws)
            if isempty(obj.s2p)
                obj.s2p = wy_wind_power_curve_data();
            end
            wp = wy_wind_speed2power(ws, obj.s2p, obj.type);
        end
    end     %% methods
end         %% classdef
