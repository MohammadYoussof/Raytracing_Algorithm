%shooting.m
% April 2022, M. Youssof!!
function [xh,zh,vh,pp,teta,time] = shooting(vpp,vps,zn,xx,xs,xr,ops)
%============================================================================================
% Some Constants
%============================================================================================
itermax = 20;
offset = abs(xs-xr); xc = 10;
%============================================================================================
% Determine Option
%============================================================================================
if ops == 1;
vh = vpp;
elseif ops == 2
vh = vps;
end
% Initial guess of the depth & time
zh = zn - 100000*eps;
t = inf*ones(size(offset));
p = inf*ones(size(offset));
%============================================================================================
% Start Raytracing
%============================================================================================
% Trial shooting
pmax = 1/min(vh);
pp = linspace(0,1/max(vh),length(xx));
sln = vh(1:length(zh)-1)*pp;
vel = vh(1:length(zh)-1)*ones(1,length(pp));
dz = abs(diff(zh))*ones(1,length(pp));
if(size(sln,1)>1)
xn = sum((dz.*sln)./sqrt(1-sln.^2));
tt = sum(dz./(vel.*sqrt(1-sln.^2)));
else
xn = (dz.*sln)./sqrt(1-sln.^2);
tt = dz./(vel.*sqrt(1-sln.^2));
end
xmax = max(xn);
%============================================================================================
% Bisection Method
%============================================================================================
% Start Bisection Method
for k=1:length(offset)
% Analyze the radius of target
n = length(xn);
xa = xn(1:n-1); xb = xn(2:n);
opt1 = xa <= offset(k) & xb > offset(k);
opt2 = xa >= offset(k) & xb < offset(k);
ind = find((opt1) | (opt2));
if (isempty(ind))
if (offset(k) >=xmax);
a = n; b = [];
else
a = []; b = 1;
end
else
a = ind; b = ind + 1;
end
x1 = xn(a); x2 = xn(b);
t1 = tt(a); t2 = tt(b);
p1 = pp(a); p2 = pp(b);
iter = 0; err= (b-a)/2;
% Minimize the error & intersect the reflector
while and(iter < itermax,abs(err) < 1)
iter = iter + 1;
xt1 = abs(offset(k) - x1); xt2 = abs(offset(k) - x2);
if and(xt1 < xc,xt1 <= xt2)
% Linear interpolation
t(k) = t1 + (offset(k)-x1)*(t2-t1)/(x2-x1);
p(k) = p1 + (offset(k)-x1)*(p2-p1)/(x2-x1);
elseif and(xt2 < xc,xt2<=xt1)
% Linear interpolation
t(k) = t2 + (offset(k)-x2)*(t1-t2)/(x1-x2);
p(k) = p2 + (offset(k)-x2)*(p1-p2)/(x1-x2);
end
% Set new ray parameter
if (isempty(a));
p2 = p1; p1 = 0;
elseif (isempty(b))
p1 = p2; p2 = pmax;
end
pnew = linspace(min([p1 p2]),max([p1 p2]),3);
% Do shooting by new ray parameter
sln = vh(1:length(zh)-1)*pnew(2);
vel = vh(1:length(zh)-1)*ones(1,length(pnew(2)));
dz = abs(diff(zh))*ones(1,length(pnew(2)));
if (size(sln,1)>1)
xtemp = sum((dz.*sln)./sqrt(1-sln.^2));
ttemp = sum(dz./(vel.*sqrt(1-sln.^2)));
else
xtemp = (dz.*sln)./sqrt(1-sln.^2);
ttemp = dz./(vel.*sqrt(1-sln.^2));
end
xnew = [x1 xtemp x2]; tnew = [t1 ttemp t2]; xmax = max(xnew);
% Analyze the radius of target
n = length(xnew);
xa = xnew(1:n-1);xb = xnew(2:n);
opt1 = xa <= offset(k) & xb > offset(k);
opt2 = xa >= offset(k) & xb < offset(k);
ind = find((opt1) | (opt2));
a = ind; b = ind + 1;
x1 = xnew(a); x2 = xnew(b);
t1 = tnew(a); t2 = tnew(b);
p1 = pnew(a); p2 = pnew(b);
err = (b - a)/2;
% Declare ray parameter
if xr > xs; pp = p; else pp = -p; end
% Compute travel time & angle
dx = real((pp.*vh.*dz)./sqrt(1-pp.*pp.*vh.*vh));
xx = xs + cumsum(dx); xh = [xs;xx];
dz = real(dx.*sqrt(1-pp.*pp.*vh.*vh))./(pp.*vh);
dt = dz./(vh.*sqrt(1-pp.*pp.*vh.*vh));
tt = cumsum(dt); time = tt(end);
teta = real(asin(pp*vh));
end % End the while loop
end % End the offset loop