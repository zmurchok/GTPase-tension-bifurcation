clear all
close all

bettercolors;

loaddata=1; %flag for plotting already computed data

if loaddata==0

  init

  global cds sys

  sys.gui.pausespecial=0;  %Pause at special points
  sys.gui.pausenever=1;    %Pause never
  sys.gui.pauseeachpoint=0; %Pause at each point

  syshandle=@SingleGTPaseTension;  %Specify system file

  SubFunHandles=feval(syshandle);  %Get function handles from system file
  RHShandle=SubFunHandles{2};      %Get function handle for ODE


  %b 0.1 1 par diagram
  bta=0;
  b=0.17;

  gam=1.5;
  GT=2;
  ell0=1;
  phi=0.75;
  Gh=0.3;
  epsilon=0.1;
  alfa=10;
  hilln=4;
  hillp=4;

  xinit=[0;1]; %Set ODE initial condition

  %Specify ODE function with ODE parameters set
  RHS_no_param=@(t,x)RHShandle(t,x,bta,b,gam,GT,ell0,phi,Gh,epsilon,alfa,hilln,hillp);

  %Set ODE integrator parameters.``
  options=odeset;
  % options=odeset(options,'RelTol',1e-12);
  options=odeset(options,'maxstep',1e-2);

  %Integrate until a steady state is found.
  [tout xout]=ode45(RHS_no_param,[0,2000],xinit,options);

  %%
  %%%%% Continuation from equilibrium %%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %Set initial condition as the endpoint of integration.  Use
  %to bootstrap the continuation.
  xinit=xout(size(xout,1),:);

  pvec=[bta,b,gam,GT,ell0,phi,Gh,epsilon,alfa,hilln,hillp]';      % Initialize parameter vector

  %active parameter
  %continue w.r.t to beta
  ap=1;
  [x0,v0]=init_EP_EP(syshandle, xinit', pvec, ap); %Initialize equilibrium
  opt=contset;
  opt=contset(opt,'MaxNumPoints',150); %Set numeber of continuation steps
  opt=contset(opt,'MaxStepsize',1e-2);
  opt=contset(opt,'Singularities',1);  %Monitor singularities
  opt=contset(opt,'Eigenvalues',1);    %Output eigenvalues
  opt = contset(opt,'InitStepsize',0.01); %Set Initial stepsize
  [x1,v1,s1,h1,f1]=cont(@equilibrium,x0,v0,opt); %GO!

  % limit cycle
  xtmp = x1(1:3,s1(2).index);
  pvec(1) = xtmp(end);
  [x0,v0] = init_H_LC(syshandle,xtmp(1:2),pvec,[1],1e-5,20,4);
  opt = contset;
  opt=contset(opt,'MaxStepsize',1e-1);%DECREASE THIS FOR BETTER ACCURACY
  opt = contset(opt,'MaxNumPoints',1500);%%INCREASE THIS FOR MORE LIMIT CYCLE
  opt = contset(opt,'Singularities',1);
  [xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);

  global lds %need this to get all the stuff for plotting
  %change the +1 after nphase to +2 for the second variable L
  %system dimension are (G,L)

  save('Fig2C.mat')
end


if loaddata==1
  load('Fig2C.mat');
end

%PLOTTING

figure
width=2.15;
set(gcf,'Units','inches','Position',[5,5,width,width*3/4])
lw = 3;

xeqcurve=x1;
minevaleq=real(f1(2,:)); %This is the last eigenvalue.  That is the one that determines stability
%
L=length(xeqcurve(1,:));
curveind=1;
lengthind=0;
maxlengthind=0;
evalstart=floor(heaviside(minevaleq(1)));
datamateq=zeros(4,L);
%
for i=1:L
    evalind=floor(heaviside(minevaleq(i)));
    if evalstart~=evalind
        curveind=curveind+1;
        i;
        evalstart=evalind;
        maxlengthind=max(lengthind,maxlengthind);
        lengthind=0;
    end
    datamateq(1,i)=xeqcurve(3,i); % This is the parameter that is varied.
    datamateq(2,i)=xeqcurve(1,i); % This is the dependent axis of the bifurcation plot.  The one you wish to plot
    datamateq(3,i)=evalind;
    datamateq(4,i)=curveind;
    lengthind=lengthind+1;
end
maxlengthind=max(maxlengthind,lengthind);
curveindeq=curveind;
for i=1:curveindeq
    index=find(datamateq(4,:)==i);
    eval(['curve' num2str(i) 'eq' '=datamateq(1:3,index);']);
end

for i=1:curveindeq
    stability=eval(['curve' num2str(i) 'eq(3,1)']);
    if stability==0
        plotsty='-';
        plotcolor=red*0.6;
        lw=2;
    else
        plotsty='-';
        plotcolor=red;
        lw=1;
    end
    plot(eval(['curve' num2str(i) 'eq(1,:)']),eval(['curve' num2str(i) 'eq(2,:)']),'color',plotcolor,'linestyle',plotsty,'Linewidth',lw)
    hold on
end

x = xlc;
s = slc;
plot(x(end,:),max(x((0:lds.tps-1)*lds.nphase+1,:)),'-','color',1-(1-green)*0.6,'linewidth',1)
plot(x(end,:),min(x((0:lds.tps-1)*lds.nphase+1,:)),'-','color',1-(1-green)*0.6,'linewidth',1)
for j = 2:size(x,2)-1 %don't do the whole curve becuase the starting point and end point shouldn't be flagged
  if ~isempty(find([s.index]==j))
    plot(x(end,j),max(x((0:lds.tps-1)*lds.nphase+1,j)),'k.')
    plot(x(end,j),min(x((0:lds.tps-1)*lds.nphase+1,j)),'k.')
  end
end

set(gca,'YColor','k','Xcolor','k','box','on','fontsize',10)
ylabel('$G$','Interpreter','latex')
xlabel('$\beta$','Interpreter','latex')
axis([0.1 0.18 0.3 1.1])
print(1,'Fig2C.eps','-depsc','-painters')
print(1,'Fig2C.png','-dpng','-r600')
