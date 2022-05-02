function [tp] = wy_wind_trans_probs(model, np, bins)

%[tp] = wy_wind_trans_probs(model, pidx0, np, bins) generates transition
%probabilities simulated using time series model for NPCC wind sites
%
% model := struct with fields:
%        ar1 = (nw_all x 1) vector of AR[1] coefficients for individual sites
%        ar1_total = scalar AR[1] coefficient for total wind
%        ols = (nw_all x 9)vector of ols estimation 
%        for [C CY1	SY1	CY2	SY2 CD1	SD1	CD2	SD2]
%        ols_total = (1 x 9) vector
%       'model_npcc.mat' => 'ar1','ar1_total','ols1','ols1_total'
%       var_wnr
%       var_wnr_total
% np    :number of periods of interest (e.g. for planning horizon)
% bins  :bin specification, supplied as either:
%        (1) number of bins (nb), or
%        (2) (1 x nb-1) vector of bin boundaries (standard deviation
%            coefficients), where initial -Inf and final +Inf are assumed,
%            but not included
%
% tp = (1 x np) cell array of transition probabilities, where
%            1st element is (1 x nb), rest are (nb x nb)
% 2022.03.27
% Wooyoung Jeon

if nargin <3
    bins = 5;
    nb=bins;
    if nargin < 2
        np = 24;
    end
end

% % test input
%  clear all;
%  load('model_npcc.mat');
%  np=24;
%  bins=5;

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


 
 rho=model.ar1_total;   %ar1 coefficient of total model
 sd_wnr = sqrt(model.var_wnr_total);    % sd of wnr of total model

 vf=[];
 sumrho=[];
var_for=[]; % var[forecasted LWIND]
for t=1:np
    sumrho(t,1) = (rho^2)^(t-1);
end

for t=1:np
    vf(t,1) = sd_wnr^2 * sum(sumrho(1:t,1));
end

sdf = sqrt(vf);
 sb=[-inf -2 -1 1 2 inf];
 sm = [];
 for i=1:bins 
     sm(i)=norminv_nst(normcdf_nst(sb(i)) + (normcdf_nst(sb(i+1)) - normcdf_nst(sb(i))) / 2) ;
 end

 tp1=[];

 midx = floor((nb+1)/2);% index for central bin

 % tp(1) : probability of being in each bin
 for j=1:nb
     tp1(midx,j,1)= normcdf_nst(sb(j+1)) - normcdf_nst(sb(j)) ; 
 end
 %sm=[-2.37322 -1.38317 0 1.38317 2.37322];
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
tp(1) = {squeeze(tp1(midx,:,1))};
for t=2:np
    tp(t) = {squeeze(tp1(:,:,t))};
end

