% logsyn.m
% April 2022, M. Youssof!!
function [pson,sson,spson,rlog,zlog] = logsyn(vp,vs,ro,z,n)
curvep = fit(z,vp,'pchipinterp'); csp = curvep.p;
curves = fit(z,vs,'pchipinterp'); css = curves.p;
curver = fit(z,ro,'pchipinterp'); csr = curver.p;
if n==0,
pson = 1.0e6./vp;
sson = 1.0e6./vs;
spson = (pson + sson)/2;
rlog = 1.0e6./ro;
zlog = z;
else
zlog = linspace(min(z),max(z),n);
plog = ppval(csp,zlog);
slog = ppval(css,zlog);
rlog = ppval(csr,zlog);
pson = 1.0e6./plog - (10)*rand(1,length(zlog));
pson = pson(:);
sson = 1.0e6./slog - (10)*rand(1,length(zlog));
sson = sson(:);
spson = (pson + sson)/2;
rlog = rlog - 0.1*rand(1,length(zlog));
rlog = rlog(:);
end;
