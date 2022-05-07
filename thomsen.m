% thomsen.m
% April 2022, M. Youssof!!
function [Rpp,incidence] = thomsen(teta1,Vp1,Vs1,rho1,teta2,Vp2,Vs2,rho2)
i = (teta1+teta2)/2;
Z1 = Vp1*rho1;
Z2 = Vp2*rho2;
Z = (Z1+Z2)/2;
dZ = Z2-Z1;
Vp = (Vp1+Vp2)/2;
dVp= Vp2-Vp1;
Vs = (Vs1+Vs2)/2;
G1 = rho1*Vs1^2;
G2 = rho2*Vs2^2;
G = (G1+G2)/2;
dG = G2-G1;
Rpp=((1/2)*(dZ/Z))+1/2*((dVp/Vp)-((2*Vs/Vp)^2*(dG/G)))*((sin(i))^2)+((1/2)*(dVp/Vp)*((sin(i))^2)*((tan(i))^2));
incidence = teta1;