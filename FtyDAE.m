function [yp] = FtyDAE( t,y )
global Int Exh QLHV SpS Runiv omega reducedCa dummy heatloss
global GaussatCA50  mfuIVCClose si EtaComb Bore N Stroke p_plenum VDisp Vc


%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Mi = [SpS.Mass];Nsp = length(Mi);
%%
pI=Int.p;TI=Int.T;
pE=Exh.p;TE=Exh.T;
%%
p=y(1);T=y(2);mi=y(3:end);
m = sum(mi);
[V,dVdt,A]=CylVolumeFie(t);
Yi = [mi/m]';
%%
Ca          = time2ang(t,omega)/(2*pi)*360;
reducedCa   = mod(Ca+360,720)-360;
CADS        = omega/(2*pi)*360;
%%
for ii=1:Nsp
    hi(ii) = HNasa(T,SpS(ii));
    ei(ii) = ENasa(T,SpS(ii));
    Cpi(ii) = CpNasa(T,SpS(ii));
    Cvi(ii) = CvNasa(T,SpS(ii));
end
h = hi*Yi';
e = ei*Yi';
Cp = Cpi*Yi';
Cv = Cvi*Yi';
gamma = Cp/Cv;
%% Intake and Exhaust thermodynamical properties
hI      = [Int.h];
eI      = [Int.e];
CpiI    = [Int.Cp];
hE      = [Exh.h];
eE      = [Exh.e];
CpiE    = [Exh.Cp];
%%
StateCyl.p          = p;
StateCyl.T          = T;
StateCyl.gamma      = gamma;
StateCyl.rho        = m/V;
StateCyl.Rg         = Cp-Cv;

StateIntake.p       = pI;
StateIntake.T       = TI;
YY                  = [Int.Y];
Ma                  = 1/(sum(YY./Mi));
Rg                  = Runiv/Ma;
StateIntake.rho     = pI/TI/Rg;
hpI                 = YY*hI';
CpI                 = YY*CpiI';
StateIntake.gamma   = CpI/(CpI-Rg);
[dmdtI]             = mdot(StateIntake,StateCyl,t,Int);


StateExhaust.p      = pE;
StateExhaust.T      = TE;
YY                  = [Exh.Y];
Ma                  = 1/(sum(YY./Mi));
Rg                  = Runiv/Ma;
hpE                 = YY*hE';       
StateExhaust.rho    = pE/TE/Rg;
CpE                 = YY*CpiE';
StateExhaust.gamma   = CpE/(CpE-Rg);
[dmdtE]             = mdot(StateExhaust,StateCyl,t,Exh);

if (abs(dmdtI) > 0)
    mfuIVCClose = mi(1);
end

if (dmdtI > 0)
    hpaI = hpI;
    YI = Int.Y;
else
    hpaI = h;
    YI = Yi;
end
if (dmdtE > 0)
    YE    = Yi;     % Same as cylinder but enthalpy is that of exhaust
    hpaE  = YE*hE';
else
    hpaE  = h;
    YE    = Yi;
    Exh.Y = YE;
end
dmidt       = [YI*dmdtI + YE*dmdtE]';
dmfuComb    = EtaComb*mfuIVCClose*pdf(GaussatCA50,reducedCa)*CADS;
dmidt_c     = si'*dmfuComb;
dmidt       = dmidt - dmidt_c;
dQcomb      = QLHV*dmfuComb;
dQcomb_real = ei*dmidt_c;

[C1, C2] = woschniconstants(reducedCa);

Twall   = 273+80;
Sp = N/60*2*Stroke;
Tr = 293;
pr = p_plenum;
Vr = Vc + VDisp;
pmotor = p_plenum*(Vc+VDisp)^gamma/(V^gamma);
w = C1*Sp + C2*(VDisp*Tr)/(pr*Vr)*(p-pmotor);
alfa    = 3.26*Bore^(-0.2)*(p/1000)^(0.8)*T^(-0.55)*w^(0.8);
dQhl        = alfa*A*(Twall-T);
heatloss(dummy) = dQhl;
dummy = dummy + 1;
%% DAE formulation
Rg = StateCyl.Rg;
yp = [dQhl-p*dVdt+hpaI*dmdtI+hpaE*dmdtE;p*V-Rg*T*m;dmidt];
end

