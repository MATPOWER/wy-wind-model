function [wsr] = wy_wind_realizations(wind_data, widx, pidx0, np)

% % wind_data := either wind_data (np x nw) with 
%        wind_speed = (nw_all x np_all) matrix of wind speeds (in m/s), corresponding
%        to nw_all specific sites, np_all periods, a particular npd, and a
%        starting dt (dt0)
%        'winddata_npcc.mat' => 'winddata': 26303 X 16
%        16 sites : ny1~ny9, ne1~ne7
%        2004, 2005, 2006 : 3 years
%        starts at 2004.01.01 1am, end at 2006.12.31 23pm
%        8760 * 3 + 24 -1 = 26303
%        wind speed data, m/s
%
%        := or model which is a struct variable with fields:
%        ar1 = (nw_all x 1) vector of AR[1] coefficients for individual sites
%        ar1_total = scalar AR[1] coefficient for total wind
%        ols = (nw_all x 9)vector of ols estimation 
%        for [C CY1	SY1	CY2	SY2 CD1	SD1	CD2	SD2]
%        ols_total = (1 x 9) vector
%       'model_npcc.mat' => 'ar1','ar1_total','ols1','ols1_total'
%       var_wnr
%       var_wnr_total
%
%   if wind_data is matrix variable, then extract realization from the data
%   if wind_data is struct variable, then generate realization
%
% widx = (nw x 1) vector of indices of wind sites of interest
% pidx0 = scalar period index of first period of horizon of interest
% wsr = (nw x np) matrix of wind speed realizations
% np    :number of periods of interest (e.g. for planning horizon)
%
% wsr : wind speed realization, (widx x np+1) : (16 x 25) 
%
% 2022.03.27
% Wooyoung Jeon

if nargin <4
    np = 24;
    if nargin <3
        pidx0 = 5112;
        if nargin <2
            widx=[1:16];
        end
    end
end

nw = length(widx);

%  clear all;
%  load('model_npcc.mat');
%  pidx0=5112;
%  np=24;
%  nw = 16;
% wind_data=model;

% if wind_data is struct variable, it is model => generate realization
if isstruct(wind_data)
model=wind_data;

 % copy lower triangular part to upper triangular part to make it symmatric
 % var-covar matrix
 var_wnr=tril(model.var_wnr,-1)'+model.var_wnr;

 % generate randomized wnr based on normal distribution using var-covar
 % matrix
 rng('default')
 gen_wnr = mvnrnd(zeros(1,nw),var_wnr,np);

 % ar(1) part
% for t=1:np
%   ar_wnr(t,:) = gen_wnr(t,:) .* model.ar1^t;
% end

for t=1:nw
    for i=1:np
        for j=1:np
            % create matrix of [e1 e1*ar1 e1*ar1^2 ... e1*ar1^23] for each
            % e_i
            temp1(i,j,t) = gen_wnr(i,t) .*model.ar1(t)^(j-1);
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
tt2=[pidx0+shift+1:1:pidx0+24+shift]';

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

for i=1:nw
    mean_logwind(:,i)= model.ols(i,1) + var_cycle * model.ols(i,2:end)'; %yfit, (np x nw)
end

 % ols + ar(1)

 realized_logwind = mean_logwind + ar_sum;

 % logwind to wind
 wsr = 10.^(realized_logwind) -1;


 % if wind_data is not a struct variable, it is actual wind speed => take
 % it from the data set
else
    wsr = wind_data(pidx0:pidx0+np,widx); 

end
