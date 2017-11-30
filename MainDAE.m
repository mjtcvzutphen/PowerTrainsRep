clear all
close all
clc
%% Add path to general functions and set Runiv
addpath('General');
global Runiv SpS QLHV
Runiv = 8.3144598;
%% Datadir just to show how you can organize output
DataDir = 'output';
%% Units, for convenience only
g=1e-3;
ms=1e-3; 
bara=1e5;
mm=1e-3;cm=1e-2;dm=0.1;
liter = dm^3;
%% Set a few global variables
global rc LCon Stroke Bore N omega Di De p_plenum dummy heatloss heatlost dt % Engine globals
LCon    = 261.6*mm;                 % connecting rod length
Stroke  = 158*mm;                   % stroke
Bore    = 130*mm;                   % bore
rc      = 17.45;                    % compression ratio
N       = 2000;                     % RPM
dummy   = 1;
Cyl.LCon = LCon;Cyl.Stroke=Stroke;Cyl.Bore=Bore;Cyl.rc=rc;
%% Simple combustion model settings (a gaussian distribution)
global GaussatCA50 mfuIVCClose si EtaComb
CA50=10;                        % CA50 (50% HR)
BDUR=20;                        % Burn Duration
GaussatCA50  = gmdistribution(CA50,BDUR,1);
EtaComb = 0.99;                 % Combustion efficiency
Comb.Shape  = GaussatCA50;
Comb.eta    = EtaComb;
Tparts = [400 475 630 415 840]; 
%% Intake and exhaust pressures
p_plenum  = 3.5*bara;               % plenum pressure
T_plenum  = 300;                    % plenum temperature
p_exhaust = p_plenum+0.1*bara;      % exhaust back-pressure
T_exhaust  = 400;                   % plenum temperature

%% Geometric and timing data of the valves
filenameValveData=fullfile('General','Valvedata.mat');
load(filenameValveData);
global Int Exh 
Di = 48.5*mm;      % diameter inlet valve
De = 43.25*mm;     % diameter exhaust valve
Int.Ca=CaI;Int.L=LI;Int.D=Di;Int.p=p_plenum;Int.T   = T_plenum;
Exh.Ca=CaE;Exh.L=LE;Exh.D=De;Exh.p=p_exhaust;Exh.T   = T_exhaust;

%% Chemistry and Thermodynamic properties
filenameThermalDataBase=fullfile('General','NasaThermDatFull.mat');
load(filenameThermalDataBase);
indexes = myfind({Sp.Name},{'Diesel','O2','N2','CO2','H2O'});   % Means that the order is now set for all vectors if you want to use SpS lateron directly
SpS     = Sp(indexes);
Nsp     = length(SpS);
Names   = {SpS.Name};
Mi      = [SpS.Mass];
Xair    = [0 0.21 0.79 0 0];                                    % Air comp
Xfuel   = [1 0 0 0 0];                                          % Fuel comp
Yair    = Xair.*Mi/(Xair*Mi');
Yfuel   = [1 0 0 0 0];
nC      = SpS(1).Elcomp(3);                                     % SpS(1) is fuel by definition (line 50)
nH      = SpS(1).Elcomp(2);
nui     = [1   nC+nH/4 0 -nC -nH/2];                            % Reaction stoichiometry (molar)
si      = nui.*Mi/Mi(1);                                        % Reaction stoichiometry (mass)
AFstoi_molar  = nui(2)+nui(2)*Xair(3)/Xair(2);                  % So-called stoichiometric air fuel ratio (fuel property for given air composition), sometimes students use this as AFstoi.
AFstoi  = si(2)+si(2)*Yair(3)/Yair(2);                          % So-called stoichiometric air fuel ratio (fuel property for given air composition)
%% Set simulation time
Ncyc    = 4;
REVS    = N/60;
omega   = REVS*2*pi;
tcyc    = (2/REVS);
t       = [0:0.1:360]./360*tcyc*Ncyc;
dt = t(2) - t(1);
%% Compute initial conditions and intake/exhaust composition
V0      = CylVolumeFie(t(1));
T0      = 273;
p0      = 3.5*bara;                                             % Typical full load set point
lambda  = 1.6;                                                  
AF      = AFstoi*lambda;                                        % Real AF ratio
EGRf    = 0.15;                                                 % EGR fraction (mass based)
fracfu  = 1/(AF+1);
fracair = 1 - fracfu;
YReactants = fracfu*Yfuel + fracair*Yair;                       % Composition vector
YProducts  = YReactants - si*YReactants(1);                     % Corresponding fully combusted setting.
for i=1:Nsp
    hi(i) = HNasa(T0,SpS(i));
end
QLHV = (hi*YReactants'-hi*YProducts')/YReactants(1);            % Classical definition of lower heating value. Just for reference not used!
Comb.QLHV    = QLHV;

Int.Y   = (1-EGRf)*YReactants+EGRf*YProducts;                   % Applying EGR setpoint
Exh.Y   = YProducts;
for ii=1:Nsp
    Int.h(ii)   = HNasa(Int.T,SpS(ii));
    Int.e(ii)   = ENasa(Int.T,SpS(ii));
    Int.Cp(ii)  = CpNasa(Int.T,SpS(ii));
    Exh.h(ii)   = HNasa(Exh.T,SpS(ii));
    Exh.e(ii)   = ENasa(Exh.T,SpS(ii));
    Exh.Cp(ii)  = CpNasa(Exh.T,SpS(ii));
end
Mave    = 1/sum([Int.Y]./Mi);
Rg      = Runiv/Mave;
mass    = p0*V0/Rg/T0;
massfu  = Int.Y(1)*mass;
Settings.N      = N;
Settings.EGR    = EGRf;
Settings.AF     = AF;
%% Set initial solution (it is an DAE problem so we must initialize)
iCase = 1;
y0(1)=p0;y0(2)=T0;y0(3:3+Nsp-1) = mass*[Int.Y];
yNames={'p','T','','','','',''};
for i=3:3+length(Names)-1
    yNames{i}=char(Names(i-2));
end
mfuIVCClose     = y0(3);
%% Solving the DAE system
tspan=t;
odopt=odeset('RelTol',1e-4,'Mass',@MassDAE,'MassSingular','yes');           % Set solver settings (it is a DAE so ...,'MassSingular','yes')
tic;
[time,y]=ode15s(@FtyDAE,tspan,y0,odopt);                                    % Take a specific solver
tel=toc;
heatlost = sum(heatloss*dt);
fprintf('Spent time %9.2f (solver %s)\n',tel,'ode15s');
%% Specify SaveName
CaseName = ['Case' num2str(iCase,'%3.3i') '.mat'];
SaveName = fullfile(DataDir,CaseName);
V = CylVolumeFie(time);
save(SaveName,'Settings','Cyl','Int','Exh','Comb','time','y','yNames','V','SpS', 'heatlost');
fprintf('Saved solution of Case %3i to %s\n',iCase,SaveName);
