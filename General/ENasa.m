function [E,EU] = ENasa(T,Sp)
global Runiv
%ENasa(T,Sp):: computes specific energy of species at Temp T
% 
%   Input: T, temperature
%          Sp, database entry
[H,~]   = HNasa(T,Sp);
E       = H-T*Runiv/Sp.Mass;
EU      = E/(Runiv/Sp.Mass);
end

