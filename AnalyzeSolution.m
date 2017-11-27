clear all;close all;clc;
%%
%% Add path to general functions and set Runiv
addpath('General');
global Runiv SpS
Runiv = 8.3144598;
%% Datadir just to show how you can organize output
DataDir = 'output';
%% Units, for convenience only
g=1e-3;
ms=1e-3; 
bara=1e5;
mm=1e-3;cm=1e-2;dm=0.1;
liter = dm^3;
%%
iCase = 1;
CaseName = ['Case' num2str(iCase,'%3.3i') '.mat'];
SaveName = fullfile(DataDir,CaseName);
load(SaveName);
fprintf('Read solution of Case %3i from %s\n',iCase,SaveName);
whos
%% Do your stuff
iSpSel = [3 4 6 7]; % Skip N2
t=time;p = y(:,1);T=y(:,2);mi=y(:,iSpSel);
RPM = Settings.N;
REVS = RPM/60;trev = 1/REVS;nREVS = (t(end)-t(1))/trev;
% Select a cycle
it = find(t > (nREVS-2)*trev & t <= nREVS*trev);
tp = t(it);
pp = p(it);
Tp = T(it);
mip = mi(it,:);
figure(1)
subplot(2,2,1)
plot(tp/ms,pp/bara);
xlabel('t [ms]');ylabel('p [bara]');
subplot(2,2,2)
plot(tp/ms,Tp);
xlabel('t [ms]');ylabel('T [K]');
subplot(2,2,[3 4])
plot(tp/ms,mip/g);
xlabel('t [ms]');ylabel('m [g]');
legend(yNames{iSpSel});
%% pV diagram
figure(2)
Vp = V(it);
subplot(1,2,1)
pl=plot(V/liter,p/bara,'--',Vp/liter,pp/bara,'r-');
set(pl(end),'LineWidth',2);
xlabel('V [l]');ylabel('p [bara]');
subplot(1,2,2)
pl=loglog(V/liter,p/bara,'--',Vp/liter,pp/bara,'r-');
set(pl(end),'LineWidth',2);
xlabel('log V [l]');ylabel('log p [bara]');
%% Computations
W   = trapz(Vp,pp); % Work, integral pdV
dummy = find(t > (nREVS-1.25)*trev); % I guess this is after IVC and before combustion. There are better ways.
index = dummy(1);
mfuel = mi(index,1);% fuel mass after intake valve close (just before combustion starts for instance)
% Plot it for checks
figure(1)
subplot(2,2,[3 4])
line(t(index)*[1 1]/ms,mfuel*[1 1]/g,'Marker','o','MarkerSize',8,'MarkerFaceColor','y');
tx=text(t(index)*[1 1]/ms,1.1*mfuel*[1 1]/g,'Selected fuel mass','Rotation',45);
QLHV = Comb.QLHV;
Qin = mfuel*QLHV;
eff = W/Qin;

%% Calcuation of MEPs
%(2.0502,3.8515)
pp_gross = pp - 3.8515*10^5;
[ja, index] = min(abs(pp_gross));
[pp_intersect,index_p] = min(abs(pp_gross));
for i=(index_p-1):size(pp_gross)
pp_gross(i) = 9999999;
end 
[pp_intersect,index_pp] = min(abs(pp_gross)); 

pp_partial = pp(index_pp:index_p);
Vp_partial = Vp(index_pp:index_p);
W_net = trapz(Vp_partial,pp_partial);
Vd = max(V) - min(V);

T = (W_net)/(2*pi);

% indicated mean effective pressure gross
IMEPg = W / Vd;
% indicated mean effective pressure net
IMEPn = W_net / Vd;
% Pump mean effective pressure
PMEP = (W-W_net) / Vd; 
% Brake mean effective pressure
BMEP = T / Vd; 
