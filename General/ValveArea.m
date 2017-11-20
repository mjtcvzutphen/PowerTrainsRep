function [A,nCyc]=ValveArea(t,CaV,LV,DV)
global omega
[ang,nCyc]=time2ang(t,omega);
Ca=ang/2/pi*360;
CaReduced = mod(Ca+360,720)-360; % Map between -360 360
% nCyc2 = floor((Ca+180)/720);
% CaL=[CaV+(nCyc-1)*720 CaV+(nCyc)*720 CaV+(nCyc+1)*720];LL=[LV LV LV];
% Lift = interp1(CaL,LL,CaReduced);
Lift = interp1(CaV,LV,CaReduced);
CurtainArea = 2*pi*DV/2*Lift;
RunnerArea  = pi*(DV/2)^2;
A = min(CurtainArea,RunnerArea);
if (isnan(Lift))
    CaL(1)/720,CaL(end)/720,Ca/720
    Ca,nCyc
    pause
end