function [Cv,CvU] = CvNasa(T,Sp)
global Runiv
%CvNasaNew(T,Sp):: computes heat capacity of species at Temp T
% 
%   Input: T, temperature
%          Sp, database entry
[~,CpU]=CpNasa(T,Sp);
CvU = CpU-1;
Cv=CvU.*Runiv/Sp.Mass;
end

