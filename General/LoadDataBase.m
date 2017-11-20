function [Spr,iOk]=LoadDataBase(Names,linfo)
global Runiv dirGeneral
Runiv = 8.3145;
addpath(dirGeneral);
%%
if (nargin==1)
linfo=false;
end
%% Load struct
load(fullfile(dirGeneral,'NasaThermDatFull'));
if (linfo)
    for i=1:length(Sp)
        fprintf('%3i :: %s\n',i,Sp(i).Name);
    end
    fprintf('Total database consist of %3i Entries,\n', length(Sp));
    uiwait(msgbox('Press ok'));
else
    fprintf('Total database consist of %3i Entries\n', length(Sp));
end
%%
[iSpecies]=myfind({Sp.Name},Names);
iOk=find(iSpecies == -1);
%% Select the subset
Spr = Sp(iSpecies);
end