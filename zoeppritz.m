% zoeppritz.m
% April 2022, M. Youssof!!
function [RC,incidence] = zoeppritz(pp,teta1,vp1,vs1,rho1,teta2,vp2,vs2,rho2,ops)
costeta1 = cos(teta1);
costeta2 = cos(teta2);
j1 = asin(pp*vs1);
j2 = asin(pp*vs2);
cosj1 = cos(j1);
cosj2 = cos(j2);
a = rho2*(1-2*vs2^2*pp^2)-rho1*(1-2*vs1^2*pp^2);
b = rho2*(1-2*vs2^2*pp^2)+2*rho1*vs1^2*pp^2;
c = rho1*(1-2*vs1^2*pp^2)+2*rho2*vs2^2*pp^2;
d = 2*(rho2*vs2^2-rho1*vs1^2);
E = b*costeta1/vp1+c*costeta2/vp2;
F = b*cosj1/vs1+c*cosj2/vs2;
G = a-d*(costeta1/vp1)*(cosj2/vs2);
H = a-d*(costeta2/vp2)*(cosj1/vs1);
D = E*F+G*H*pp^2;
if ops==1
Rpp = ((b*costeta1/vp1-c*costeta2/vp2)*F-(a+d*(costeta1/vp1)*(cosj2/vs2))*H*pp^2)/D;
RC = Rpp;
elseif ops==2
Rps = (-2*(costeta1/vp1)*(a*b+c*d*(costeta2/vp2)*(cosj2/vs2)*(pp*vp1)))/(vs1*D);
RC = Rps;
end
incidence = teta1;