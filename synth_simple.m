% model.m
% april 2022 , M. Youssof
%==========================================================================================
% MAIN PROGRAM
%==========================================================================================
clc; clear all; close all;
format short g
%============================================================================================
% Define Geometry
%============================================================================================
fprintf('---> Defining the interfaces ...\n');
% Input geometry
xmin = 0; xmax = 3000;
zmin = 0; zmax = 300;
% Make strata layer
zlayer = [0;40; 80;180;300];%;450;630;800;1050;1200;1400;1620;1820;2000];
nlayer = length(zlayer);
layer = 1:1:nlayer;
thick = abs(diff(zlayer));
x = [xmin xmax]; z = [zlayer zlayer];
dg = 10;
xx = xmin:dg:xmax; nx = length(xx);
zz = repmat(zlayer,1,nx);
fprintf(' Geometry has been defined...[OK]\n');
%============================================================================================
% Source-Receiver Groups
%============================================================================================
fprintf('---> Setting source-receiver configuration ...\n');
% Receiver Interval
dr = 20;
% Mix Configuration
% Source
% xs = [55 1350];
% zs = [0 690];
% ns = length(xs);
%
% % Receiver
% xr = [200 500 850 1600 1900];
% zr = [300 0 220 150 0];
% nr = length(xr);
% --------------------------------------
% Front Configuration
% Source
% xs = xmin;
% zs = 0;
% ns = length(xs);
%
% % Receiver
% xr = xmin+dr:dr:xmax;
% zr = zeros(length(xr));
% nr = length(xr);
% --------------------------------------
% End Configuration
% Source
% xs = xmax;
% zs = 0;
% ns = length(xs);
%
% % Receiver
% xr = xmax-dr:-dr:xmin;
% zr = zeros(length(xr));
% nr = length(xr);
% --------------------------------------
% Split Spread Configuration
% Source
% xs = fix((xmax-xmin)/2);
% zs = 0;
% ns = length(xs);
%
% % Receiver
% xr = [fix((xmax-xmin)/2)-dr:-dr:xmin fix((xmax-xmin)/2)+dr:dr:xmax];
% zr = zeros(length(xr));
% nr = length(xr);
% --------------------------------------
% Off-End Configuration
% Source
% xs = [xmin xmax];
% zs = [0 0];
% ns = length(xs);
%
% % Receiver
% xr = xmin+dr:dr:xmax-dr;
% zr = zeros(length(xr));
% nr = length(xr);
% Front Spread Configuration
% Source
xs = [xmin fix((xmax-xmin)/2) xmax];
zs = [0 0 0];
ns = length(xs);
% Receiver
xr = [xmin+dr:dr:fix((xmax-xmin)/2)-dr fix((xmax-xmin)/2)+dr:dr:xmax-dr];
zr = zeros(length(xr));
nr = length(xr);
% VSP Configuration
% Source
% xs = 1950;
% zs = 0;
% ns = length(xs);
%
% % Receiver
% zr = 0:100:1600;
% xr = 50.*ones(length(zr));
% nr = length(xr);
nray = ns*nr;
fprintf(' Source-Receiver Groups have been setted...[OK]\n');
%============================================================================================
% Create Elastic Parameter
%============================================================================================
fprintf('---> Creating elastic parameter ...\n');
% Create synthetic Vp, Vs, and Density
vlayer = [1800;1300;2000;2200;2500];%;2400;2700;3150;2950;3440;3750;4000;4350;4600]; % Velocity
vlayer = vlayer(1:length(zlayer));
vel = [vlayer vlayer];
vp = vlayer; % P wave velocity
vs = (vp-1360)./1.16; % S wave velocity based on Castagna's rule
ro = 0.31.*(vp).^0.25; % Density based on Gardner's rule
pois = (vs.^2-0.5*vp.^2)./(vs.^2-vp.^2); % Poisson ratio
fprintf(' No.Layer Depth(m) Vp(m/s) Vs(m/s) Density(kg/m^3) Poisson Ratio \n');
disp([layer' zlayer vp vs ro pois]);
fprintf(' Elastic parameters have been created...[OK]\n');
figure;
set(gcf,'color','white');
% Plot P wave velocity
subplot(1,3,1);
stairs(vp,zlayer,'LineWidth',2.5,'Color',[0.07843 0.1686 0.549]);
ylabel('Depth (m)','FontWeight','bold','Color','black');
title('V_{p} (m/s)','FontWeight','bold');
set(gca,'XMinorGrid','on','YMinorGrid','on','YDir','reverse','YColor',[0.04314 0.5176 0.7804],...
'XAxisLocation','top','XColor',[0.04314 0.5176 0.7804],'MinorGridLineStyle','-','FontWeight',...
'demi','FontAngle','italic');
% Plot S wave velocity
subplot(1,3,2);
stairs(vs,zlayer,'LineWidth',2.5,'Color',[1 0 0]);
title('V_{s} (m/s)','FontWeight','bold');
set(gca,'XMinorGrid','on','YMinorGrid','on','YDir','reverse','YColor',[0.04314 0.5176 0.7804],...
'XAxisLocation','top','XColor',[0.04314 0.5176 0.7804],'MinorGridLineStyle','-','FontWeight',...
'demi','FontAngle','italic');
% Plot Density
subplot(1,3,3);
stairs(ro,zlayer,'LineWidth',2.5,'Color',[0.07059 0.6392 0.07059]);
title('Density (kg/m^3)','FontWeight','bold');
set(gca,'XMinorGrid','on','YMinorGrid','on','YDir','reverse','YColor',[0.04314 0.5176 0.7804],...
'XAxisLocation','top','XColor',[0.04314 0.5176 0.7804],'MinorGridLineStyle','-','FontWeight',...
'demi','FontAngle','italic');
% Plot Geology Model
figure;
set(gcf,'color','white');
%mohammad color = load('mycmap.txt');

color = colormap(copper); colorbar ('horz');

pcolor(xx,zz,repmat(vlayer,1,nx)); shading flat; hold on
colormap(color); colorbar ('horz');
axis([0 xmax zmin-0.03*zmax zmax])
set(gca,'YDir','reverse','XaxisLocation','bottom',....
'Ytick',zlayer,'FontWeight','demi','PlotBoxAspectRatioMode','Manual',...
'PlotBoxAspectRatio',[2.4 1.2 1],'Position',[0.04 0.30 0.90 0.60]);
% Plot Source-Receiver Group
plot(xs,zs,'r*','markersize',12); hold on
plot(xr,zr,'sk','markersize',4,'markerfacecolor','c'); hold on
xlabel('Velocity (m/s)','FontWeight','bold','Color','black');
ylabel('Depth (m)','FontWeight','bold','Color','black');
%============================================================================================
% Run Ray Tracing
%============================================================================================
fprintf('---> Starting ray tracing ...\n');
wat = waitbar(0,'Raytracing is being processed, please wait...');
%============================================================================================
xoff = [];
% Loop over for number of souce
tic
for i=1:ns
% Loop over for number of receiver
for j=1:nr
%====================================================================================
% Loop over for number of layer
for k=1:nlayer
%=====================================================================================
% Declare reflection boundary
if and(zr(j) < zlayer(k),zs(i) < zlayer(k))
zm = zz(k,:); zf = min(zm)- 100000*eps;
% Downgoing path
d = find(zlayer > zs(i));
if(isempty(d)); sdown = length(zlayer); else sdown = d(1)-1; end
d = find(zlayer > zf);
if(isempty(d)); edown = length(zlayer); else edown = d(1)-1; end
zd = [zs(i);zlayer(sdown+1:edown)]; nd = length(zd);
% Upgoing path
u = find(zlayer > zr(j));
if(isempty(u)); sup = length(zlayer); else sup = u(1)-1; end
u = find(zlayer > zf);
if(isempty(u)); eup = length(zlayer); else eup = u(1)-1; end
zu = [zr(j);zlayer(sup+1:eup+1)]; nu = length(zu);
zn = [zd;(flipud(zu))]; nrefl = length(zn)-1;
%=================================================================================
% Declare elastic parameter
% Downgoing elastic parameter
vpd = [vp(sdown:edown);vp(edown)];
vsd = [vs(sdown:edown);vs(edown)];
rod = [ro(sdown:edown);ro(edown)];
% Upgoing elastic parameter
vpu = [vp(sup:eup);vp(eup)];
vsu = [vs(sup:eup);vs(eup)];
rou = [ro(sup:eup);ro(eup)];
% Combine model elastic parameter
vpp = [vpd(1:end-1);flipud(vpu(1:end-1))];
vss = [vsd(1:end-1);flipud(vsu(1:end-1))];
vps = [vpd(1:end-1);flipud(vsu(1:end-1))];
rho = [rod(1:end-1);flipud(rou(1:end-1))];
%=================================================================================
% Start Raytracing (P-P, S-S, or P-S mode)
ops = 1; % ops=1 for PP mode; ops=2 for PS mode
[xh,zh,vh,pp,teta,time] = shooting(vpp,vps,zn,xx,xs(i),xr(j),ops);
theta = abs(teta); twt(k,j,i) = time;
% Plot Ray
if ops == 1
plot(xh,zh,'k-');
title('Seismic Raytracing (P-P mode)','FontWeight','bold');
elseif ops == 2
xd = xh(1:nd+1); xu = xh(nd+1:end);
zd = zh(1:nd+1); zu = zh(nd+1:end);
plot(xd,zd,'k-',xu,zu,'r-');
title('Seismic Raytracing (P-S mode)','FontWeight','bold');
end
%=================================================================================
% Compute Reflection Coefficient (Downgoing-Upgoing)
for c=1:nrefl-1
% Reflection Coefficient of Zoeppritz Approximation
[rc1,teta] = zoeppritz(pp,theta(c),vpp(c),vss(c),rho(c),theta(c+1),...
vpp(c+1),vss(c+1),rho(c+1),ops);
rcz(c,j,i) = rc1;
% Reflection Coefficient of Shuey Approximation
[A,B,rc2,teta] = shuey(theta(c),vpp(c),vss(c),rho(c),theta(c+1),...
vpp(c+1),vss(c+1),rho(c+1),ops);
rcs(c,j,i) = rc2; AA(c,j,i) = A; BB(c,j,i) = B;
% Reflection Coefficient of Thomsen Approximation
[rc3,teta] = thomsen(theta(c),vpp(c),vss(c),rho(c),theta(c+1),...
vpp(c+1),vss(c+1),rho(c+1));
rct(c,j,i) = rc3;
angle(c,j,i) = teta.*(180/pi);
end
%====================================================================================
end
end % for horizon/reflector end
xoff = [xoff xr(j)];
waitbar(j/nray,wat)
end % for receiver end
% Save Data
save(['time_shot',num2str(i),'.mat'],'twt');
save(['reflz_shot',num2str(i),'.mat'],'rcz');
save(['refls_shot',num2str(i),'.mat'],'rcs');
save(['reflt_shot',num2str(i),'.mat'],'rct');
save(['teta_shot',num2str(i),'.mat'],'angle');
save(['intercept_shot',num2str(i),'.mat'],'AA');
save(['gradient_shot',num2str(i),'.mat'],'BB');
waitbar(i/ns,wat)
end % for sources end
%=============================================================================================
close(wat);
toc
fprintf(' Ray tracing has succesfully finished...[OK]\n');
xx = repmat(xr,nlayer,1);
% Plot Traveltime
figure;
for i=1:ns
tm = load(['time_shot',num2str(i),'.mat']); tt = tm.twt;
plot(xx(2:nlayer,:),tt(2:nlayer,:,i),'b.'); hold on
xlabel('Horizontal Position (m)','FontWeight','bold','FontAngle','normal','Color','black');
ylabel('Time (s)','FontWeight','bold','FontAngle','normal','Color','black');
title('Traveltime','FontWeight','bold');
set(gca,'Ydir','reverse');
set(gcf,'color','white');
end
% Plot Reflection Coefficient
figure;
for i=1:ns
reflz = load(['reflz_shot',num2str(i),'.mat']); reflz = reflz.rcz;
refls = load(['refls_shot',num2str(i),'.mat']); refls = refls.rcs;
reflt = load(['reflt_shot',num2str(i),'.mat']); reflt = reflt.rct;
theta = load(['teta_shot',num2str(i),'.mat']); angle = theta.angle;
%Reflection Coefficient of Zoeppritz
subplot(1,3,1)
plot(abs(angle(1:nlayer,:,i)),reflz(1:nlayer,:,i),'r.'); hold on
xlabel('Incidence Angles (degree)','FontWeight','bold','Color','black');
ylabel('Reflection Coefficient','FontWeight','bold','Color','black');
grid on; title('Rpp Zoeppritz','FontWeight','bold','Color','black')
set(gca,'YColor',[0.04314 0.5176 0.7804],'XColor',[0.04314 0.5176 0.7804]); hold on
%Reflection Coefficient of Shuey
subplot(1,3,2)
plot(abs(angle(1:nlayer,:,i)),refls(1:nlayer,:,i),'g.'); hold on
xlabel('Incidence Angles (degree)','FontWeight','bold','Color','black');
grid on; title('Rpp Shuey','FontWeight','bold','Color','black')
set(gca,'YColor',[0.04314 0.5176 0.7804],'XColor',[0.04314 0.5176 0.7804]); hold on
%Reflection Coefficient of Thomsen
subplot(1,3,3)
plot(abs(angle(1:nlayer,:,i)),reflt(1:nlayer,:,i),'b.'); hold on
xlabel('Incidence Angles (degree)','FontWeight','bold','Color','black');
grid on; title('Rpp Thomsen','FontWeight','bold','Color','black')
set(gca,'YColor',[0.04314 0.5176 0.7804],'XColor',[0.04314 0.5176 0.7804]); hold on
set(gcf,'color','white');
end
for i=1:ns
refls = load(['refls_shot',num2str(i),'.mat']); refls = refls.rcs;
theta = load(['teta_shot',num2str(i),'.mat']); angle = theta.angle;
Rt = load(['intercept_shot',num2str(i),'.mat']); Rt = Rt.AA;
Gt = load(['gradient_shot',num2str(i),'.mat']); Gt = Gt.BB;
Ro(1:nlayer,:,i) = Rt(1:nlayer,:,i);
Go(1:nlayer,:,i) = Gt(1:nlayer,:,i);
Rc(1:nlayer,:,i) = refls(1:nlayer,:,i);
inc(1:nlayer,:,i) = angle(1:nlayer,:,i);
end
Rp = reshape(Ro,nr*ns*nlayer,1);
G = reshape(Go,nr*ns*nlayer,1);
R = reshape(Rc,nr*ns*nlayer,1);
teta = reshape(inc,nr*ns*nlayer,1);
% Plot Attributes
figure;
% Rp-G Cross Plot
subplot(1,2,1)
plot(Rp,G,'r.'); hold on
%i dont know 
% plot(Rl,G,'m-'); hold on
xlabel('R','FontWeight','bold','Color','black');
ylabel('G','FontWeight','bold','Color','black');
grid on; title('R-G Cross Plot','FontWeight','bold','Color','black')
set(gca,'YColor',[0.04314 0.5176 0.7804],'XColor',[0.04314 0.5176 0.7804]); hold on
% R-sin^2(teta) Cross Plot
subplot(1,2,2)
plot((sin(teta).^2),R,'g.'); hold on
xlabel('sin^2(teta)','FontWeight','bold','Color','black');
ylabel('R(teta)','FontWeight','bold','Color','black');
grid on; title('R-sin^2(teta) Cross Plot','FontWeight','bold','Color','black')
set(gca,'YColor',[0.04314 0.5176 0.7804],'XColor',[0.04314 0.5176 0.7804]); hold on
set(gcf,'color','white');
%============================================================================================
% AVO Modelling
%============================================================================================
fprintf('---> Starting AVO Modelling ...\n');
% Make wavelet ricker (f = 2 Hz)
dt = 0.004; f = 20;
[w,tw] = ricker(dt,f);
% Create zeros matrix for spike's location
tmax = max(tt(:));
tr = 0:dt:tmax; nt = length(tr);
t = tt(2:nlayer,:);
% Take reflectivity into spike location
spikes = zeros(nt,nray);
rcz = real(reflz(2:nlayer,:));
for k=1:nlayer-1
for j=1:nray
ir(k,j) = fix(t(k,j)/dt+0.1)+1;
spikes(ir(k,j),j) = spikes(ir(k,j),j) + rcz(k,j);
end
end


%%%%%
%%%%%%
%%%%
%%%%% Convolve spikes with wavelet ricker
for j=1:nray
seisz(:,j) = convz(spikes(:,j),w);
%seisz = conv2(spikes(:,j),w);
end
 ampz=reshape(seisz,length(tr),nr,ns);
% ampz=(permute(seisz,[3 2 1]));
for i=1:ns
ampz_shot(:,:,i) = ampz(:,:,i);
save(['ampz_shot',num2str(i),'.mat'],'ampz_shot');
end
% CMP Stacked
Stacked = zeros(nt,nr);
for i=1:ns
files = load(['ampz_shot',num2str(i),'.mat']);
files = files.ampz_shot;
Stacked = Stacked + files(:,:,i);
end
% Save Seismic Data
save(['seis','.mat'],'Stacked');
% Plot AVO
figure;
wig(xr,tr,Stacked,'black'); hold on
xlabel('Offset (m)','FontWeight','bold','Color','black');
ylabel('Time (s)','FontWeight','bold','Color','black');
title('Synthetic Seismogram','FontWeight','bold','Color','black');
axis tight
set(gca,'YColor',[0.04314 0.5176 0.7804],'XColor',[0.04314 0.5176 0.7804]); hold on
set(gcf,'color','white');
fprintf(' AVO Modelling has been computed...[OK]\n');
%============================================================================================
% Seismogram Synthetic
%============================================================================================
% Create Smooth Synthetic Log From Elastic Data Set
[pson,sson,spson,rlog,zlog] = logsyn(vp,vs,ro,zlayer,500);
vplog = 1.0e6./pson; vslog = 1.0e6./sson;
rholog = 0.31.*(vplog).^0.25; vpvs = vplog./vslog;
poisson = (vslog.^2-0.5*vplog.^2)./(vslog.^2-vplog.^2);
% Plot Synthetic Data Log
figure;
set(gcf,'color','white');
% Plot P Sonic Log
subplot(1,3,1);
stairs(pson,zlog,'LineWidth',1,'Color',[0.07843 0.1686 0.549]);
xlabel('P sonic (US/m)','FontWeight','bold','Color','black');
ylabel('Depth (m)','FontWeight','bold','Color','black');
set(gca,'XMinorGrid','on','YMinorGrid','on','YDir','reverse','XDir','reverse','YColor',...[0.04314 0.5176 0.7804],'XAxisLocation','top','XColor',...
[0.04314 0.5176 0.7804],'MinorGridLineStyle','-','FontWeight','demi','FontAngle','italic');
% Plot S Sonic Log
subplot(1,3,2);
stairs(sson,zlog,'LineWidth',1,'Color',[1 0 0]);
xlabel('S sonic (US/m)','FontWeight','bold','Color','black');
set(gca,'XMinorGrid','on','YMinorGrid','on','YDir','reverse','XDir','reverse','YColor',[0.04314 0.5176 0.7804],...
'XAxisLocation','top','XColor',[0.04314 0.5176 0.7804],'MinorGridLineStyle','-','FontWeight',...
'demi','FontAngle','italic');
% Plot Density Log
subplot(1,3,3);
stairs(rlog,zlog,'LineWidth',1,'Color',[0.07059 0.6392 0.07059]);
xlabel('Density (kg/m^3)','FontWeight','bold','Color','black');
set(gca,'XMinorGrid','on','YMinorGrid','on','YDir','reverse','YColor',[0.04314 0.5176 0.7804],...
'XAxisLocation','top','XColor',[0.04314 0.5176 0.7804],'MinorGridLineStyle','-','FontWeight',...
'demi','FontAngle','italic');
% Plot Rock Properties
figure;
% Cross Plot P and S wave velocity
subplot(2,2,1);
plot(vplog,vslog,'.','Color',[0.07843 0.1686 0.549]);
xlabel('Vp (m/s)','FontWeight','bold','Color','black');
ylabel('Vs (m/s)','FontWeight','bold','Color','black');
title('Cross Plot of Vp and Vs','FontWeight','bold','Color','black');
set(gca,'XMinorGrid','on','YMinorGrid','on','YColor',[0.04314 0.5176 0.7804],...
'XAxisLocation','bottom','XColor',[0.04314 0.5176 0.7804],'MinorGridLineStyle',...
'-','FontWeight','demi','FontAngle','italic');
% Plot S Sonic Log
subplot(2,2,2);
plot(vplog,vpvs,'.','Color',[1 0 0]);
xlabel('Vp (m/s)','FontWeight','bold','Color','black');
ylabel('Vp/Vs ratio','FontWeight','bold','Color','black');
title('Cross Plot of Vp and Vp/Vs Ratio','FontWeight','bold','Color','black');

% set(gca,'XMinorGrid','on','YMinorGrid','on','YColor',[0.04314 0.5176 0.7804],...
%'XAxisLocation','bottom','XColor',[0.04314 0.5176 0.7804],...
%'MinorGridLineStyle','','FontWeight',...
%'demi','FontAngle','italic');

% Plot Density Log
subplot(2,2,3);
plot(rholog,vpvs,'.','Color',[0.07059 0.6392 0.07059]);
xlabel('Density (kg/m^3)','FontWeight','bold','Color','black');
ylabel('Vp/Vs ratio','FontWeight','bold','Color','black');
title('Cross Plot of Vp and Vp/Vs Ratio','FontWeight','bold','Color','black');
set(gca,'XMinorGrid','on','YMinorGrid','on','YColor',[0.04314 0.5176 0.7804],...
'XAxisLocation','bottom','XColor',[0.04314 0.5176 0.7804],'MinorGridLineStyle',...
'-','FontWeight','demi','FontAngle','italic');
subplot(2,2,4);
plot(poisson,vpvs,'m.');
xlabel('Poisson Ratio','FontWeight','bold','Color','black');
ylabel('Vp/Vs ratio','FontWeight','bold','Color','black');
title('Cross Plot of Vp and Vp/Vs Ratio','FontWeight','bold','Color','black');
set(gca,'XMinorGrid','on','YMinorGrid','on','YColor',[0.04314 0.5176 0.7804],...
'XAxisLocation','bottom','XColor',[0.04314 0.5176 0.7804],'MinorGridLineStyle',...
'-','FontWeight','demi','FontAngle','italic');
set(gcf,'color','white');
% Integrate Depth and Velocity To Get Two-Way Time
tstart = 0; tracenum = 10;
[tz,zt] = sonic2tz(spson,zlog,-100,tstart);
tzobj = [zt tz];
vins = vplog; dens = rholog; z = zlog;
