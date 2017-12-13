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
it = find(t > (nREVS-2)*trev & t <= nREVS*trev);    % taking the last two revolutions (last cycle)
tp = t(it);                                         % using only these values for time, pressure, Temperature etc. 
pp = p(it);
Tp = T(it);
mip = mi(it,:);
Vp = V(it);
% Plot a cycle's t/p, t/T and t/m diagrams
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
subplot(1,2,1)
pl=plot(V/liter,p/bara,'--',Vp/liter,pp/bara,'r-');
set(pl(end),'LineWidth',2);
xlabel('V [l]');ylabel('p [bara]');
subplot(1,2,2)
pl=loglog(V/liter,p/bara,'--',Vp/liter,pp/bara,'r-');
set(pl(end),'LineWidth',2);
xlabel('log V [l]');ylabel('log p [bara]');
%% Computations
W   = trapz(Vp,pp);                  % total work, integral pdV
dummy = find(t > (nREVS-1.25)*trev); % I guess this is after IVC and before combustion. There are better ways.
index = dummy(1);
mfuel = mi(index,1); % fuel mass after intake valve close (just before combustion starts for instance)
% Plot it for checks
figure(1)
subplot(2,2,[3 4])
line(t(index)*[1 1]/ms,mfuel*[1 1]/g,'Marker','o','MarkerSize',8,'MarkerFaceColor','y');
tx=text(t(index)*[1 1]/ms,1.1*mfuel*[1 1]/g,'Selected fuel mass','Rotation',45);
QLHV = Comb.QLHV;
Qin = mfuel*QLHV;
eff = W/Qin;

%% Calcuation of MEPs
% self-determined cutting coordinates: [2.0465,3.8597]

% Determining index numbers in which the intersection point certainly is located
index_legit_V = find(V > 1.95*10^-3);   % getting the indexes that surely have the intersection point
pp_gross = pp - 3.8597*10^5;            % moving the point of intersect to p=0 in order to later use the min() function to find the index of the intersect

% cutting off the parts that surely don't contain the intersection point
pp_gross_cut=[];                        % initializing the vector
for i = 1:length(index_legit_V);
    tempindex = index_legit_V(i);       % using the index values from the index_legit_V vector instead of value i
    pp_gross_cut(tempindex) = pp_gross(tempindex);  % copying all the "pp_gross values inside the index domain" to the new ``cut-vector''
end
% changing the 0 values to NaN values
for i = 1:length(pp_gross_cut);
    if pp_gross_cut(i) == 0
        pp_gross_cut(i) = NaN;
    else
    end
end

% determining first intersection index number
[pp_intersect1,index_pp1] = min(abs(pp_gross_cut));   % finding the index of p in that intersection point in the cut vector

% new for-loop by Menno Z
for i=(index_pp1-100):(index_pp1+100);
pp_gross_cut(i) = 9999999;
end 

% determining second intersection index number
[pp_intersect2,index_pp2] = min(abs(pp_gross_cut-pp_intersect1)); 

pp_partial = pp(index_pp1:index_pp2);
Vp_partial = Vp(index_pp1:index_pp2);
W_net = trapz(Vp_partial,pp_partial);
Vd = max(V) - min(V);           % determining the volume inside the
figure(2)
subplot(1,2,1)
hold on
plot(Vp_partial/liter,pp_partial/bara,'b-')
T = (W_net)/(2*pi);             % determining torque

% indicated mean effective pressure gross
IMEPg = W / Vd;
% indicated mean effective pressure net
IMEPn = W_net / Vd;
% Pump mean effective pressure
PMEP = (W-W_net) / Vd; 
% Brake mean effective pressure
BMEP = T / Vd; 
