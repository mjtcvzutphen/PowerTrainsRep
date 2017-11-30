function [dY] = SingleWiebe(CA,a,n,dCA,CAI)
% Computes single wiebe (vibe) function for given parameter set a,n,dCA,CAI
% at crankangle CA
%   [dY] = SingleWiebe(CA,a,n,dCA,CAI) 
%   Input:  CA,a,n,dCA,CAI
%           arg = ((CA-CAI)/dCA)
%           dY = a*n/dCA * arg^(n-1)*exp(-a*arg^n)
    ARG = (CA-CAI)./dCA;
    ARG = max(ARG,0);   % Because of possible fractional powers of n this must be imposed. Otherwise imaginary numbers.
    EXPARG = -a*ARG.^(n);
    PF = ARG.^(n-1);
    EXPF=exp(EXPARG);
    dY = a*n/dCA*(PF.*EXPF);
end

