function [IntOut, CumOut]=Integr(x,y)
% Integr [IntOut, CumOut]=Integr(x,y)
%      Integrates Int y dx
%         Intout: Integral
%         Cumout(i): Int_x(1)^x(i) y dx 
Nx=length(x);
Sz=size(x);
dx=(x(2)-x(1))/2;
S=y(1)*dx;
CS(1)=S;
for i=2:Nx-1
    dx=(x(i)-x(i-1));
    S=S+(y(i)+y(i-1))/2*dx;
    CS(i)=S;
end;
dx=(x(Nx)-x(Nx-1))/2;
S=S+y(Nx)*dx;
CS(Nx)=S;
Sz2=size(CS);
if (Sz2(1) ~= Sz(1)) 
    CS=CS';
end;
IntOut=S;CumOut=CS;

