function [PET,VPD] = penmont_model_q(Ta,Rnet,ea,press)

% [PET] = penmont_model_q(Ta,Rnet,e,press)
%
% Modified to take data from a model
%
% This function will use the Penman-Monteith method to calculate reference
% (potential) evapotranspiration from DAILY data.  Radiation requirements
% are calculated in situ, based on elevation, latitude, tmax, and tmin.
% 
% PET  = potential evapotranspiration (mm day-1) 
%
% Ta        = temperature, degrees C
% e         = Vapor pressure, Pascals
% Rnet      = surface net radiation, W/m2
% press     = surface pressure, Pa
% leap      = 1 if leap year, 0 if not
%
% Written by Benjamin I. Cook
% 
% Based on: Xu and Singh (2002), from Water Resources Management
%           Lawrence (2005), BAMS
%
%           FAO Document:
%           http://www.fao.org/docrep/X0490E/x0490e00.htm#Contents
%
% For Tetens, both above and below zero:
%       http://cires.colorado.edu/~voemel/vp.html

% Constants
B1=243.04;  % Magnus formulation
A1=17.625;  % Magnus formulation

coast=0;
leap=0;
u=1;

if coast==1
    distcoast=0.19;
elseif coast==0
    distcoast=0.16;
end

% Calculate the latent heat of vaporization (MJ kg-1)
lambda_lv=2.501-(2.361e-3).*Ta;

% Convert Pressure to kPa
press=press./1000;

% Saturation vapor pressure (kPa)
es=0.611.*exp((17.27.*Ta)./(Ta+237.3));   

% Actual vapor pressure (kPa)
%ea=0.611.*exp((17.27.*Td)./(Td+237.3));   
ea=ea./1000;

% Dewpoint Temperature (degrees C)
%Td=((Ta+273.15).*(1-(((Ta+273.15).*log(rh./100))./(2.501e6./461.5))).^(-1))-273.15;

% Relative Humidity
rh=ea./es;

% Slope of the vapor pressure curve (kPa C-1)
%delta_vpc=(4098.*es.*Ta)./((Ta+237.3).^2);
delta_vpc=(4098.*es)./((Ta+237.3).^2);

% Psychometric constant (kPa C-1)
psych_const=0.00163.*(press./lambda_lv);

% Net radiation
Rnet=Rnet./11.6;  % convert W/m2 to MJ/m2/d

% Soil Heat Flux
%    soil heat capacity*daily temp change*effective soil depth
%gflux=0.7.*tdiff.*0.15;
%    at daily time steps, we can ignore ground heat flux
gflux=0;

VPD=es-ea;

% Potential Evapotranspiration (mm)
PET=(0.408.*delta_vpc.*(Rnet-gflux)+psych_const.*(900./(Ta+273)).*u.*(es-ea))./(delta_vpc+psych_const.*(1+0.34.*u));









