% shuey.m
% April 2022, M. Youssof!!
function [A,B,RC,incidence] = shuey(teta1,Vp1,Vs1,rho1,teta2,Vp2,Vs2,rho2,ops)
i=(teta1+teta2)/2;
Z1 = Vp1*rho1;
Z2 = Vp2*rho2;
Z = (Z1+Z2)/2;
dZ = Z2-Z1;
Vp = (Vp1+Vp2)/2;
dVp = Vp2-Vp1;
Vs = (Vs1+Vs2)/2;
dVs = Vs2-Vs1;
Rb = dVs/(2*Vs);
gam = (Vs1+Vs2)/(Vp1+Vp2);
rho = (rho1+rho2)/2;
Rro = (rho2-rho1)/(2*rho);
po1 = (((Vp1/Vs1)^2)-2)/(2*(((Vp1/Vs1)^2)-1));
po2 = (((Vp2/Vs2)^2)-2)/(2*(((Vp2/Vs2)^2)-1));
po = (po1+po2)/2;
dpo = po2-po1;
A = (1/2)*(dZ/Z);
B = -2*(1-(po/(1-po)))*A-(1/2*((1-3*po)/(1-po))*(dVp/Vp))+(dpo/(1-po)^2);
C = 1/2*(dVp/Vp);
if ops==1
Rpp = A+B*(sin(i))^2+C*(sin(i))^2*(tan(i))^2;
RC = Rpp;
elseif ops==2
Rps = -gam*((Rro+2*gam(2*Rb+Rro))*sin(i));
RC = Rps;
end
incidence = teta1;