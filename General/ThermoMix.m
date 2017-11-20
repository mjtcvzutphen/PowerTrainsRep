function [Cp,Cv,H,E] = ThermoMix(Y,T,Sp)
global Runiv
%
% Y = squeeze(Yf);
Ma=0;
for i=1:length(Sp);
    Hi(i) = HNasaNew(T,Sp(i));
    Cpi(i)= CpNasaNew(T,Sp(i));
    Ma    = Ma+Y(i)/Sp(i).Mass;
end
Ma = 1/Ma;
H  = Y*Hi';  
E  = H-Runiv/Ma*T;
Cp = Y*Cpi';
Cv = Cp- Runiv/Ma;
end

