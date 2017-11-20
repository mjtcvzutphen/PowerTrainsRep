function [phim] = Bernoulli(State1,State2)
%Bernoulli's function
%   Input: State1 and State2 containing the states at both sides of the
%   orifice
%   Output: isentropic mass flow rate per unit area.
%           NO Discharge Coeff used here
Ru = 8.314;
dpstate=(State2.p - State1.p);dwidth=State1.p/1000;
sgn=-tanh(dpstate/dwidth);
if (State1.p > State2.p)
%     sgn=+1;
    pu=State1.p;Tu=State1.T;rhou=State1.rho;gamma=State1.gamma;
    pt=State2.p;
else
%     sgn=-1;
    pu=State2.p;Tu=State2.T;rhou=State2.rho;gamma=State2.gamma;
    pt=State1.p;    
end
TT0=2/(gamma+1);
prchoke=(TT0)^(gamma/(gamma-1));
practual=pt/pu;
pr=max(practual,prchoke);
flag = (pr==practual);
rho=rhou*(pr)^(1/gamma);
Usq=2*gamma/(gamma-1)*(pu/rhou)*(1-pr^((gamma-1)/gamma));
phim=rho*sqrt(Usq)*sgn;
% a0sq=gamma*pu/rhou;a0=sqrt(a0sq);
% phimch=rhou*a0*sqrt((TT0)^((gamma+1)/(gamma-1)));


