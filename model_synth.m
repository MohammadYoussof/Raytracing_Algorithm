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

% Front Spread Configuration
% Source
xs = [xmin fix((xmax-xmin)/2) xmax];
zs = [0 0 0];
ns = length(xs);
% Receiver
xr = [xmin+dr:dr:fix((xmax-xmin)/2)-dr fix((xmax-xmin)/2)+dr:dr:xmax-dr];
zr = zeros(length(xr));
nr = length(xr);

nray = ns*nr;
fprintf(' Source-Receiver Groups have been setted...[OK]\n');
%============================================================================================
% Create Elastic Parameter
%============================================================================================
fprintf('---> Creating elastic parameter ...\n');
% Create synthetic Vp, Vs, and Density
vlayer = [1500;1600;2000;2200;2500];%;2400;2700;3150;2950;3440;3750;4000;4350;4600]; % Velocity
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

shuttleVideo = VideoReader('shuttle.avi');
workingDir = tempname;
mkdir(workingDir)
mkdir(workingDir,'images')

% Loop over for number of sources
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

 %imwrite(x,sprintf('%d.pdf',j))

% 
%  while hasFrame(shuttleVideo)
%    img = readFrame(shuttleVideo);
%    filename = [sprintf('%03d',j) '.jpg'];
%    fullname = fullfile(workingDir,'images',filename);
%    imwrite(img,fullname)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
%    j = j+1;
%  end
% 
%  imageNames = dir(fullfile(workingDir,'images','*.jpg'));
% imageNames = {imageNames.name}';
% 
% 
% outputVideo = VideoWriter(fullfile(workingDir,'shuttle_out.avi'));
% outputVideo.FrameRate = shuttleVideo.FrameRate;
% open(outputVideo)
end % for receiver end
% Save Data
save(['time_shot',num2str(i),'.mat'],'twt');
save(['reflz_shot',num2str(i),'.mat'],'rcz');
save(['refls_shot',num2str(i),'.mat'],'rcs');
save(['reflt_shot',num2str(i),'.mat'],'rct');
save(['teta_shot',num2str(i),'.mat'],'angle');
save(['intercept_shot',num2str(i),'.mat'],'AA');
save(['gradient_shot',num2str(i),'.mat'],'BB');
% 
% 
% for i=1:10
    imwrite(x,sprintf('%d.jpg',i))
% end;
waitbar(i/ns,wat)
end % for sources end

close(wat);
toc
fprintf(' Ray tracing has succesfully finished...[OK]\n');
xx = repmat(xr,nlayer,1);