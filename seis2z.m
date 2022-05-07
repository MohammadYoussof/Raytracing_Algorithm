% seis2z.m
% April 2022, M. Youssof!!
function [seisz,z,dz,tz,zt]=seis2z(seis,t,vp,zv,vs,ops,zmin,zmax)
% converts PP or PS surface seismic data from time to depth.
%
% seis : PP or PS data in time.
% t : time axis of PP or PS data.
% vp : P-wave velocity log.
% zv : depth axis for velocity log.
% vs : S-wave velocity log.
% ops : conversion type, ops=1 for PP data and ops=2 for PS data
% zmin : minimum desired depth.
% zmax : maximum desired depth.
% seisz : the output PP or PS data.
% z : the output depth-axis.
if ops==1; % For PP data
vs = vp;
end;
v = vs; % Get the same sample for both PP and PS data.
dt = t(2)-t(1); % Time interval
[nt,nx] = size(seis); % nt is the number of time steps.
pson=(1.e06./vp); % P-sonic
sson=(1.e06./vs); % S-sonic
spson=(pson + sson)/2; % Fake PS-sonic
% Computes an approximate 2-way time-depth curve from a sonic log for use with depth conversion.
tstart=0;
[tz,zt] = sonic2tz(spson,zv,-10,tstart); % Get time-depth curve.
dz = min(v)*dt/2; % Sample depth
nz = floor((zmax-zmin)/dz)+1; % Number of depth sample
z = (0:nz-1)*dz + 0; % Depth axis
z = z(:); % Vector for depth
t2 = interp1(zt,tz,z);
nt2 = length(t2);
seisz = zeros(nt2,nx); % Data in depth, nx = number of offsets
for k=1:nx,
seisz( : ,k) = sinci(seis(: ,k),t',t2); % Interpolate the amplitude at each depths sample
end
