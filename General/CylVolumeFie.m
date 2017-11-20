function [V,dVdt,A] = CylVolumeFie(t)
global  rc LCon Stroke Bore N omega % Engine globals
VDisp   = pi*(Bore/2)^2*Stroke;
Vc      = VDisp/(rc-1);             % clearance volume (m^3)
Rc      = Stroke/2;                 % Radius crank
r       = Rc/LCon;
radian  = time2ang(t,omega);

PistonArea  = pi*(Bore/2)^2;

x   = Rc*cos(radian)+LCon*sqrt(1-(r*sin(radian)).^2);             % Displacement w/r TDC
V   = Vc+(LCon+Rc-x)*PistonArea;                                  % Volume 
Vp  = -Rc*(N/30)*pi*sin(radian).*(1+(r*cos(radian))./(sqrt(1-(r*sin(radian)).^2)));                               
                                                                            % Piston speed
dVdt = -Vp.*pi*(Bore/2)^2;
height = V/PistonArea;
LinerArea   = height*pi*Bore;
A    = LinerArea+2*PistonArea;

end