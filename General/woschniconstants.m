function [C1, C2] = woschniconstants(crankangle)
% This function defines C1 and C2 as dictated in Heywood p.680
if (crankangle >= -360 && crankangle < -180) % intake
    C1 = 6.18;
    C2 = 0;
elseif crankangle >= -180 && crankangle < 0 % compression
    C1 = 2.28;
    C2 = 0;
elseif crankangle >= 0 && crankangle < 180 % expansion/combustion
    C1 = 2.28;
    C2 = 3.24e-3;
elseif crankangle >= 180 && crankangle < 360 % exhaust
    C1 = 6.18;
    C2 = 0;
end
end