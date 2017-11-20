function [M] = MassDAE( t,y )
global SpS
% Mass matrix belonging to MainDAE
%   Determines M of M dydt= Fty;
%   Assumes : y = [p T m1 ... m5]; Fty = [EnConv; Eq of State;dm1dt;...;dm5dt];
p=y(1);T=y(2);mi=y(3:end); %m=sum(mi);
V=CylVolumeFie(t);
Nsp = length(SpS);
Cvi     = mi; % Will be overwritten
ei      = mi; % idem
for ii=1:Nsp
    Cvi(ii) = CvNasa(T,SpS(ii));
    ei(ii) = ENasa(T,SpS(ii));
end
mCv = Cvi'*mi;           % m * Cv
%%
M = [0 mCv ei';...
     0   0  0 0 0 0 0;...
     0   0  1 0 0 0 0;...
     0   0  0 1 0 0 0;...
     0   0  0 0 1 0 0;...
     0   0  0 0 0 1 0;...
     0   0  0 0 0 0 1;...
     ];
end

