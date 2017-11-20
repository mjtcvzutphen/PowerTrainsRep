function [mdot] = mdot(S1,S2,t,V)
% Computes massflowrate using ValveArea and Bernoulli's equation
%   S1 and S2 are the states at point 1 and 2
%   V contains valve data (V.Ca,V.L,V.D)
mdot = Bernoulli(S1,S2)*ValveArea(t,V.Ca,V.L,V.D);
end

