function [wp] = wy_wind_speed2power(ws, s2p)

% ws = (np x nb x nw) 3D matrix of wind speeds (ex, 24 x 5 x 16)
%
% s2p = i) ((ws_max+1) x 1) vector (to be used for all sites of interest)
%       of wind power output percentage for matching wind speed [0:30]
%       ii) scalar to select a type of wind power curve from
%       'WindPowerCurveIEC.txt'
%        'WindPowerCurveIEC.txt' has 5 power curve,  (0 m/s ~ 30 m/s)
%        1: IEC1, 2:IEC2, 3:IEC3, 4:Offshore, 5:Multi-turbine
%        default is 5:Multi-turbine
%       used as lookup table for converting wind speed to wind power, as
%       fraction of installed capacity, where 1st and last cols correspond
%       to wind speeds of 0 m/s and ws_max m/s, respectively, 
%
% wp = (nw x np x nb) 3D matrix of available wind power outputs as
%             fraction of installed capacity
% 2022.04.01
% Wooyoung Jeon

if nargin <2
    s2p=5;
end

% test input
% s2p=5;
% ws=wsf;

[np,nb,nw]=size(ws); % dimension from ws

powercurvefile='WindPowerCurveIEC.txt'; % call windpowercurve data
powercurve=load(powercurvefile);

%if s2p is scalar, use a powercurve from the file and select based on s2p value(1~5) 
%if s2p is vector, use the entered powercurve
if isscalar(s2p)
    idx_ws =  powercurve(:,1);
    pc = powercurve(:,s2p+1);
else
    pc = s2p;
end

%% Variables from input files

% To increase precision of interpolation, increase wind speed step from 1
% to 0.01
e=(numel(pc))*100+1;
pgeneration=zeros(e,1);

for i=1:(e)
    pgeneration(i)=interp1(idx_ws,pc,(i-1)*0.01);
    if isnan(pgeneration(i))==1
        pgeneration(i)=0;
    end
end

% interpolate ws with respect to idx_ws_100 and pgeneration
idx_ws_100 = [0:0.01:max(idx_ws)+1];
wp=interp1(idx_ws_100,pgeneration,ws);
