function [ang,nCyc] = time2ang(t,om)
%Relation between time and ang
ang     = -2*pi+om*t; % choose t=0 to coincide with -2*pi which is TDC
tCyc    = 2/(om/(2*pi));
nCyc    = floor(t/tCyc);
end

