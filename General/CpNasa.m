function [Cp,CpU] = CpNasa(T,Sp)
global Runiv
%CpNasaNew(T,Sp):: computes heat capacity of species at Temp T
% 
%   Input: T, temperature
%          Sp, database entry
if (isnan(Sp.Ts))
    Tl=T;
    a=Sp.Pol(1,:);
    CpU=a(1)+a(2).*Tl+a(3).*Tl.^2+a(4).*Tl.^3+a(5).*Tl.^4;        % Formula 5.4 of lecture notes
else
    ilow = (T <= Sp.Ts);
    Tl=T(ilow);
    a=Sp.Pol(1,:);
    CpU(ilow)=a(1)+a(2).*Tl+a(3).*Tl.^2+a(4).*Tl.^3+a(5).*Tl.^4;        % Formula 5.4 of lecture notes
    
    ihigh = (T > Sp.Ts);
    Tl=T(ihigh);
    a=Sp.Pol(2,:);
    CpU(ihigh)=a(1)+a(2).*Tl+a(3).*Tl.^2+a(4).*Tl.^3+a(5).*Tl.^4;        % Formula 5.4 of lecture notes
end
Cp=CpU.*Runiv/Sp.Mass;
end

