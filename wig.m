% wig.m
% April 2022, M. Youssof!!
function wig(x,t,Data,style,dmax,showmax,plImage,imageax)
cla;
showmax_def = 500;
style_def = 'wiggle';
if nargin==1,
Data = x;
t = 1:1:size(Data,1);
x = 1:1:size(Data,2);
dmax=max(Data(:));
style=style_def;
showmax=showmax_def;
end
if nargin==2,
Data = x; dmax = t;
t = 1:1:size(Data,1);
x = 1:1:size(Data,2);
style=style_def;
showmax=showmax_def;
end
if nargin==3,
style = style_def;
dmax = max(abs(Data(:)));
showmax = showmax_def;
end
if nargin==4,
dmax=max(abs(Data(:)));
showmax=showmax_def;
end
if nargin==5,
showmax=showmax_def;
end
if nargin<7
plImage=0;
end
if isempty(dmax),
dmax=max(abs(Data(:)));
end
if isempty(showmax),
showmax=100;
end
if nargin==7,
imageax=[-1 1].*dmax;
end
if plImage==1,
imagesc(x,t,Data);
if (length(imageax)==1)
imageax=[-1 1].*abs(imageax);
end
caxis(imageax);
hold on
end
if (showmax>=0)
if length(x)>1
dx = x(2)-x(1);
else
dx = max(Data)*15;
end
ntraces = length(x);
d = ntraces/showmax;
if d<=1; d=1; end
d=round(d);
dmax = dmax/d;
for i=1:d:length(x)
xt=dx*Data(:,i)'./dmax;
if (strmatch('black',style)==1)
xt1 = xt; xt1(find(xt1>0))=0;
f1=fill(x(i)+[xt,fliplr(xt1)],[t,fliplr(t)],[0 0 0]); hold on
set(f1,'LineWidth',0.0001)
set(f1,'EdgeColor',[0 0 0])
plot(xt+x(i),t,'k-','linewidth',.05); hold on
elseif (strmatch('red',style)==1)
xt2 = xt; xt2(find(xt2>0))=0;
f2=fill(x(i)+[xt,fliplr(xt2)],[t,fliplr(t)],[1 0 0]); hold on
set(f2,'LineWidth',0.0001)
set(f2,'EdgeColor',[1 0 0])
plot(xt+x(i),t,'r-','linewidth',.05); hold on
elseif (strmatch('green',style)==1)
xt3 = xt; xt3(find(xt3>0))=0;
f3=fill(x(i)+[xt,fliplr(xt3)],[t,fliplr(t)],[0 1 0]); hold on
set(f3,'LineWidth',0.0001)
set(f3,'EdgeColor',[0 1 0])
plot(xt+x(i),t,'g-','linewidth',.05); hold on
elseif (strmatch('blue',style)==1)
xt4 = xt; xt4(find(xt4>0))=0;
f4=fill(x(i)+[xt,fliplr(xt4)],[t,fliplr(t)],[0 0 1]); hold on
set(f4,'LineWidth',0.0001)
set(f4,'EdgeColor',[0 0 1])
plot(xt+x(i),t,'b-','linewidth',.05);hold on
else
plot(xt+x(i),t,'k-','linewidth',.05); hold on
end
if i==1, hold on;end
end
end
if length(x)>1
axis([min(x)-1.2*x(1) max(x)+1.2*x(1) min(t) max(t)])
else
axis([min(x)-1.2*min(x) max(x)+1.2*max(x) min(t) max(t)])
end
set(gca,'Ydir','reverse')