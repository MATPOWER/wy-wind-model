function [wsf] = wy_wind_forecasts(model, widx, pidx0, ws0, np, bins)

% % wind_data := struct with fields:
%        wind_speed = (nw_all x np_all) matrix of wind speeds (in m/s), corresponding
%        to nw_all specific sites, np_all periods, a particular npd, and a
%        starting dt (dt0)
%        wind_cap = (nw_all X 1) :matrix of wind capacity (in MW)
%        'winddata_npcc.mat' => 'winddata': 26303 X 16, 'windcap': 1 X 16
%        16 sites : ny1~ny9, ne1~ne7
%        2004, 2005, 2006 : 3 years
%        starts at 2004.01.01 1am, end at 2006.12.31 23pm
%        8760 * 3 + 24 -1 = 26303
%        wind speed data, m/s
% widx = (nw x 1) vector of indices of wind sites of interest
% pidx0 = scalar period index of first period of horizon of interest
% ws0 = (nw x 1), initial wind speed
% np    :number of periods of interest (e.g. for planning horizon)
% bins  :bin specification, supplied as either:
%        (1) number of bins (nb), or
%        (2) (1 x nb-1) vector of bin boundaries (standard deviation
%            coefficients), where initial -Inf and final +Inf are assumed,
%            but not included
%
% wsf : wind speed forecasts, (np x nb x nw) : (ex, 25 x 5 x 16)
%       first hour : realization, later hours : forecasts
%
% 2022.03.27
% Wooyoung Jeon

if nargin <6
    bins = 5;
    nb=bins;
    if nargin < 5
        np = 24;
        if nargin <4
            pidx0 = 5112;
            if nargin <3
                widx=[1:16];
            end
        end
    end
end


if isscalar(bins)
    nb=bins;
    if nb==1
        bin_bound = [-inf inf];
    elseif nb==2
        bin_bound = [-inf 0 inf];
    elseif nb==3
        bin_bound = [-inf -1 1 inf];
    elseif nb==4
        bin_bound = [-inf -1 0 1 inf];
    elseif nb==5        
        bin_bound = [-inf -2 -1 1 2 inf];
    elseif nb==6        
        bin_bound = [-inf -2 -1 0 1 2 inf];
    elseif nb==7        
        bin_bound = [-inf -3 -2 -1  1 2 3 inf];
    elseif nb==8        
        bin_bound = [-inf -3 -2 -1 0 1 2 3 inf];
    elseif nb==9        
        bin_bound = [-inf -4 -3 -2 -1 1 2 3 4 inf];
    elseif nb==10        
        bin_bound = [-inf -4 -3 -2 -1 0 1 2 3 4 inf];
    else
        bin_bound = [-inf -2 -1 1 2 inf]; % default value if bins<0 or bins>10
    end
else
    nb=length(bins)+1;
    bin_bound = [-inf bins inf];
end

probcum = normcdf_nst(bin_bound);   % cdf for each sd
probref = probcum(1:end-1) + diff(normcdf_nst(bin_bound))/2; % find mean prob location for each bin in normal dist
bin_mean = norminv_nst(probref); % find mean location for each bin in normal dist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. initial setup for dataset and variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nw = length(widx);
nb=bins;

%winddata = log_winddata;

% apply estimates for total wind
coef_cycle = model.ols(:,2:end); % coefficient for cycles, (8x1)
coef_mean = model.ols(:,1);    % coefficient for mean(constant), (1x1)

rho = model.ar1;
sd_wnr= sqrt(model.var_wnr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. inputs from econometric model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cycle information defined
% PERIODS
PY1 = 8766;
PY2 = PY1 / 2;
PD1 = 24;
PD2 = PD1 / 2;

% adjustment for calender cycle input, push 1 hour to the next
% IMPORTANT: cycle hour, 1hour shifted as estimation is done this way
% hour 0 to hour 24, 25hours, 
% hour 0 needed for ar(1) process. at(t-1) is needed
shift = 1;
tt2=[pidx0+shift:1:pidx0+24+shift]';

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. generates forecasted bin, ff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% yy:  yact : realized LWIND
% ww: yfit : ols fitted LWIND => MEAN
% uu: yres : ols residual of LWIND
% yres_fit : ar(1) fitted LWIND
% ee: wnr : ar(1) white noise residual of LWIND
% ff: yfor : forecasted LWIND = yfit+yres_fit = ww + yres_fit


for i=1:nw

yINIT = ws0(i);
ww= coef_mean(widx(i)) + var_cycle * coef_cycle(widx(i),:)'; %yfit
uINIT = yINIT - ww(1);


sumrho=[];
var_for=[]; % var[forecasted LWIND]
for t=1:np
    sumrho(t,1) = (rho(widx(i))^2)^(t-1);
end

for t=1:np
    var_for(t,1) = sd_wnr(widx(i), widx(i))^2 * sum(sumrho(1:t,1));
end

sd_for = sqrt(var_for); % sqrt(var[forecasted LWIND])


% forecasted LWIND
ff0 = zeros(np+1,1);
ff0(1) = yINIT;
for t=1:np
    ff0(t+1) = ww(t+1) - rho(widx(i))*ww(t) + rho(widx(i))*ff0(t);
end
ff = ff0(2:end);

bm = ff0(2:end) + sd_for * bin_mean;

% uuu : AR residual for bin mean i in hour t-1
% bf : forecast for bin mean i in hour t

bf0=[];
bf=[];

uuu0 = zeros(np+1,nb); % t starting at 0, 25 hours

midx = floor((nb+1)/2);% index for central bin

uuu0(1,midx) = uINIT;
for t=2:np+1
    uuu0(t,:) = bm(t-1,:) - ff0(t,1);
end

uuu = uuu0(2:end,:);

for t=2:np
    bf0(t-1,:) = ff(t,1) + rho(widx(i))*uuu0(t,:);
end

bf_temp = ones(1,nb) * ff(1);
bf = [bf_temp;bf0];

bf_all(:,:,i) = bf; 
wsf(:,:,i) = bf;
end

%wsf = (10.^wsf_log)-1; % convert log(wind+1) tp wind speed