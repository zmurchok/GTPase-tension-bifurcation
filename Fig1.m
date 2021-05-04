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

  bta=0.135;
  b=0.146;
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
  opt=contset(opt,'MaxNumPoints',80); %Set numeber of continuation steps
  opt=contset(opt,'MaxStepsize',1e-2);
  opt=contset(opt,'Singularities',1);  %Monitor singularities
  opt=contset(opt,'Eigenvalues',1);    %Output eigenvalues
  opt = contset(opt,'InitStepsize',0.01); %Set Initial stepsize
  [x1,v1,s1,h1,f1]=cont(@equilibrium,x0,v0,opt); %GO!

  figure
  hold on
  cpl(x1,v1,s1,[3,1])%quick plot

  pause

  % limit cycle
  xtmp = x1(1:3,s1(2).index);
  pvec(1) = xtmp(end);
  [x0,v0] = init_H_LC(syshandle,xtmp(1:2),pvec,[1],1e-5,20,4);
  opt = contset;
  opt=contset(opt,'MaxStepsize',1e-2);
  opt = contset(opt,'MaxNumPoints',2000);%%INCREASE THIS FOR MORE LIMIT CYCLE
  opt = contset(opt,'Singularities',1);
  [xlc,vlc,slc,hlc,flc]=cont(@limitcycle,x0,v0,opt);


  x = xlc;
  s = slc;
  global lds %need this to get all the stuff for plotting
  %chagne the +1 after nphase to +2 for the second variable L
  %system dimension are (G,L)
  plot(x(end,:),max(x((0:lds.tps-1)*lds.nphase+1,:)),'linewidth',1)
  plot(x(end,:),min(x((0:lds.tps-1)*lds.nphase+1,:)),'linewidth',1)
  for j = 2:size(x,2)-1 %don't do the whole curve becuase the starting point and end point shouldn't be flagged
    if ~isempty(find([s.index]==j))
      plot(x(end,j),max(x((0:lds.tps-1)*lds.nphase+1,j)),'r.')
      plot(x(end,j),min(x((0:lds.tps-1)*lds.nphase+1,j)),'r.')
    end
  end

  pause
  print(1,'1param.png','-dpng','-r600')
  close



  %ANOTHER OPTION FOR THE HOMOCLINIC CONTINUATION IS TO START FROM A NEARLY HOMOCLINIC ORBIT AND THEN USE THAT! THIS MAY BE EASIER TO DO.... USING SINGULARITIES ON HERE WILL THROW AN ERROR ABOUT LINE 160 of HOMOCLINIC ...
  % aps = [1,2];
  % p(ap) = xlc(end,end);
  % T = xlc(end-1,end)/2;
  % [x0,v0] = init_LC_Hom(syshandle,xlc(:,end),slc(:,end),pvec,aps,50,4,[0,1,1],T,0.001,0.001);
  % opt = contset;
  % opt = contset(opt,'MaxStepsize',0.001);
  % opt = contset(opt,'Singularities',1); %THIS BREAKS IT FOR SOME REASON!
  % opt = contset(opt,'MaxNumPoints',200);
  % [xh,vh,sh,hh,fh] = cont(@homoclinic,x0,v0,opt);
  % figure
  % hompts2 = xh(end-5:end-4,:);
  % plot(hompts2(1,:),hompts2(2,:),'color',red,'linewidth',1)
  % %






  %%
  %%%%%%%%%%%%%%%%
  % codim 2 bif
  %b vs beta
  aps = [1,2];
  %


  pause
  %get hopf point, continue hopt point
  xtmp = x1(1:end-1,s1(2).index);
  pvec(ap) = x1(end,s1(2).index);
  [x0,v0] = init_H_H(syshandle,xtmp,pvec,aps);
  opt=contset;
  opt=contset(opt,'MaxNumPoints',1600);
  opt=contset(opt,'Singularities',1);
  opt=contset(opt,'MaxStepsize',1e-8);
  opt=contset(opt,'backward',1);
  [x2b,v2b,s2b,h2b,f2b]=cont(@hopf,x0,v0,opt);

  xtmp = x1(1:end-1,s1(2).index);
  pvec(ap) = x1(end,s1(2).index);
  [x0,v0] = init_H_H(syshandle,xtmp,pvec,aps)
  opt=contset;
  opt=contset(opt,'MaxNumPoints',1600); %Set numeber of continuation steps
  opt=contset(opt,'Singularities',1);
  opt=contset(opt,'MaxStepsize',1e-8);
  opt=contset(opt,'backward',0);
  [x2,v2,s2,h2,f2]=cont(@hopf,x0,v0,opt);

  figure()
  hold on
  cpl(x2,v2,s2,[3,4])
  cpl(x2b,v2b,s2b,[3,4])
  % xlim([0 0.5])
  % ylim([0 0.3])
  % grid


  %BELOW MAKES THE LPC CURVE (GREEN IN THE FIGURE... IT'S  VERY VERY VERY VERY SLOW, so I would just recommend loading the data I already computed and plotting from there (SEE BELOW)....)
  % pause
  % %
  % [x0,v0]=init_GH_LPC(syshandle,x2b,pvec,s2b(2),aps,30,4,0.001);
  % opt = contset;
  % opt = contset(opt,'MaxNumPoints',12000);
  % opt=contset(opt,'Singularities',1);
  % % opt=contset(opt,'MaxStepsize',0.01);
  % opt=contset(opt,'MaxStepsize',0.1);
  % [x5,v5,s5,h5,f5]=cont(@limitpointcycle,x0,v0,opt); %takes 2185 seconds
  % [x5,v5,s5,h5,f5] = cont(x5,v5,s5,h5,f5,cds); %extend
  % plot(x5(end-1,:),x5(end,:),'color',green,'linewidth',1)
  % for j = 2:size(x5,2)-1 %don't do the whole curve becuase the starting point and end point shouldn't be flagged
  %   if ~isempty(find([s5.index]==j))
  %     plot(x5(end-1,j),x5(end,j),'r-',color)
  %   end
  % end

  %
  % [x0,v0]=init_CPC_LPC(syshandle,xlc,slc(2),aps,20,4);
  % opt=contset;
  % opt=contset(opt,'MaxNumPoints',5000);
  % opt=contset(opt,'Singularities',1);
  % opt = contset(opt,'MaxStepsize',1e-2);
  % opt = contset(opt,'MaxNewtonIters',3);
  % opt = contset(opt,'MaxCorrIters',10);
  % opt = contset(opt,'MaxTestIters',10);
  % opt = contset(opt,'VarTolerance',1e-6);
  % opt = contset(opt,'FunTolerance',1e-5);
  % opt = contset(opt,'Adapt',3);
  % opt = contset(opt,'Backward',1);
  % [xlpc,vlpc,slpc,hlpc,flpc] = cont(@limitpointcycle,x0,v0,opt);
  %
  %last two rows of xlpc are the parameters
  % plot(xlpc(end-1,:),xlpc(end,:))
  % for j = 2:size(xlpc,2)-1 %don't do the whole curve becuase the starting point and end point shouldn't be flagged
  %   if ~isempty(find([slpc.index]==j))
  %     plot(xlpc(end-1,j),xlpc(end,j),'r.')
  %   end
  % end




  pause

  disp('SADDLE NODE 1')

  %get LP1
  xtmp = x1(1:end-1,s1(3).index);
  pvec(ap) = x1(end,s1(3).index);
  [x0,v0] = init_LP_LP(syshandle,xtmp,pvec,aps)
  opt = contset(opt,'MaxNumPoints',1000);
  %backward
  opt=contset(opt,'backward',1);
  [x3b,v3b,s3b,h3b,f3b]=cont(@limitpoint,x0,v0,opt);
  %forward
  opt=contset(opt,'backward',0);
  [x3,v3,s3,h3,f3]=cont(@limitpoint,x0,v0,opt);
  %plots
  cpl(x3,v3,s3,[3,4])
  cpl(x3b,v3b,s3b,[3,4])


  disp('SADDLE NODE 2')

  %get LP2
  xtmp = x1(1:end-1,s1(5).index);
  pvec(ap) = x1(end,s1(5).index);
  [x0,v0] = init_LP_LP(syshandle,xtmp,pvec,aps)
  %backward
  opt=contset(opt,'backward',1);
  [x4b,v4b,s4b,h4b,f4b]=cont(@limitpoint,x0,v0,opt);
  %forward
  opt=contset(opt,'backward',0);
  [x4,v4,s4,h4,f4]=cont(@limitpoint,x0,v0,opt);
  %plots
  cpl(x4,v4,s4,[3,4])
  cpl(x4b,v4b,s4b,[3,4])
  %
  %
  disp('HOMOCLINIC')



  %%BT point to start homoclinic
  %%% THIS IS INCREDIBLY SENSITIVE TO NUMERICS?
  %s2(3)
  %GLbetab %OLD
  xtmp = x2(1:2,s2(2).index)
  pvec(aps) = x2(3:4,s2(2).index)
  %initial
  [x0,v0] = init_BT_Hom(syshandle,xtmp,s2(2),pvec,aps,50,4,0.00001,0.0001,[0 1 1]);
  opt=contset;
  % opt=contset(opt,'Singularities',1);
  opt=contset(opt,'MaxStepsize',0.01);
  opt = contset(opt,'MaxNumPoints',500); %use 3000 to get full curve
  [xlc2,vlc2,slc2,hlc2,flc2] = cont(@homoclinic,x0,v0,opt);
  hompts = xlc2(end-5:end-4,:);
  plot(hompts(1,:),hompts(2,:),'color',purple,'linewidth',1)

  %or can continue using the Hom_Hom....
  pvec(ap) = xlc2(end-5,end-4,:); %update parameter vector !!!ap only here not aps!!!
  T = slc2(end).data.T
  [x0,v0] = init_Hom_Hom(syshandle,xlc2,vlc2,slc2(end),pvec,aps,50,4,[0 1 1],T,0.00001,0.0001);
  opt=contset;
  % opt=contset(opt,'Singularities',1); %This breaks sometimes...
  opt=contset(opt,'MaxNumPoints',70);
  opt=contset(opt,'InitStepsize',0.01);
  opt=contset(opt,'MaxStepsize',0.01);
  [xh,vh,sh,hh,fh] = cont(@homoclinic,x0,v0,opt);

  hompts2 = xh(end-5:end-4,:);
  plot(hompts2(1,:),hompts2(2,:),'color',green,'linewidth',2)

end

if loaddata==1
  load('Fig1.mat');
end


figure()
hold on
%SN
sn = plot(x3(end-1,:),x3(end,:),'color',red*0.6,'linewidth',1)
plot(x3b(end-1,:),x3b(end,:),'color',red*0.6,'linewidth',1)
plot(x4(end-1,:),x4(end,:),'color',red*0.6,'linewidth',1)
snic = plot([x4(end-1,1:515),0.165097],[x4(end,1:515),0.140843],'color',red*0.5+yellow*0.5,'linewidth',1)
plot(x4b(end-1,:),x4b(end,:),'color',red*0.6,'linewidth',1)

%NS
% ns = plot(x2(end-2,s2(2).index:end),x2(end-1,s2(2).index:end),'-','color',red,'linewidth',1)

%Hopf
hopf = plot(x2(end-2,1:s2(2).index),x2(end-1,1:s2(2).index),'color',green*0.6,'linewidth',1)
plot(x2b(end-2,:),x2b(end-1,:),'color',green*0.6,'linewidth',1)

%LPC
lpc = plot(x5(end-1,:),x5(end,:),'color',green*1.2,'linewidth',1)

%Hom
hompts = xlc2(end-5:end-4,:);
hom = plot(hompts(1,:),hompts(2,:),'-','color',yellow,'linewidth',1)


%HOM-NS DOT (approx)
homns=[0.16462912684,0.1419595605];
%SNIC-H DOT (approx)
snich=[0.165097,0.140843];
plot(homns(1),homns(2),'k.')
plot(snich(1),snich(2),'k.')

set(gca,'YColor','k','Xcolor','k','box','on','fontsize',10)
ylabel('$b$','Interpreter','latex')
xlabel('$\beta$','Interpreter','latex')

%loop over curves and add bifurction points
%lpc
for j = 2:size(x5,2)-1 %loop over curve
  if ~isempty(find([s5.index]==j)) %flag cusp points
    plot(x5(end-1,j),x5(end,j),'k.') %label with black dots
  end
end
%hopfs
for j = 2:size(x2,2)-1 %loop over curve
  if ~isempty(find([s2.index]==j)) %flag points
    plot(x2(end-2,j),x2(end-1,j),'k.') %label with black dots
  end
end
for j = 2:size(x2b,2)-1 %loop over curve
  if ~isempty(find([s2b.index]==j)) %flag points
    plot(x2b(end-2,j),x2b(end-1,j),'k.') %label with black dots
  end
end
%SNs
for j = 2:size(x3,2)-1 %loop over curve
  if ~isempty(find([s3.index]==j)) %flag points
    plot(x3(end-1,j),x3(end,j),'k.') %label with black dots
  end
end
for j = 2:size(x3b,2)-1 %loop over curve
  if ~isempty(find([s3b.index]==j)) %flag points
    plot(x3b(end-1,j),x3b(end,j),'k.') %label with black dots
  end
end
%SNs
for j = 2:size(x4,2)-1 %loop over curve
  if ~isempty(find([s4.index]==j)) %flag points
    plot(x4(end-1,j),x4(end,j),'k.') %label with black dots
  end
end
for j = 2:size(x4b,2)-1 %loop over curve
  if ~isempty(find([s4b.index]==j)) %flag points
    plot(x4b(end-1,j),x4b(end,j),'k.') %label with black dots
  end
end

set(gcf,'Units','inches','Position',[5,5,3.25,3.25*3/4])
%zoom in
% axis([0.1520    0.1700    0.1340    0.1484])
% print(1,'2par_zoom_in.png','-dpng','-r300')
% print(1,'2par_zoom_in.eps','-depsc','-painters')

axis([0.152    0.1700    0.1340    0.148])
legend('hide')
print(1,'Fig1B.png','-dpng','-r600')
print(1,'Fig1B.eps','-depsc','-painters')


yline(0.12,'k:')
yline(0.145,'k:')
yline(0.17,'k:')

% legend([sn,hopf,hom,ns,lpc,snic],{'SN','H','HOMC','NS','SNP','SNIC'},'NumColumns',2,'fontsize',8)
legend([sn,hopf,hom,lpc,snic],{'SN','H','HOMC','SNP','SNIC'},'NumColumns',1,'fontsize',8)


axis([0.08    0.22        0.1    0.25])
print(1,'Fig1A.png','-dpng','-r600')
print(1,'Fig1A.eps','-depsc','-painters')



% %zoom way out
% axis([0.0520,    0.8895,   -0.1943,    0.2245])
% print(1,'2par_zoom_out.png','-dpng','-r300')
% print(1,'2par_zoom_out.eps','-depsc','-painters')









%%% A DECENT 1D PLOTTING SCRIPT... from BILL HOLMES + other code that I'm not using is down here...

%%%%%% NICE Plotting script.  x=continuation info  f=eigenvalues %%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% width=3;
% height=3;
% x0 = 5;
% y0 = 5;
% fontsize = 12;
% f = figure('Units','inches','Position',[x0 y0 width height],'PaperPositionMode','auto');
% figure
% lw = 1;
% colordef(gcf,'black')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xeqcurve=x1;
% minevaleq=real(f1(2,:)); %This is the last eigenvalue.  That is the one that determines stability
% %
% L=length(xeqcurve(1,:));
%
% curveind=1;
% lengthind=0;
% maxlengthind=0;
% evalstart=floor(heaviside(minevaleq(1)));
% datamateq=zeros(4,L);
% %
% for i=1:L
%     evalind=floor(heaviside(minevaleq(i)));
%     if evalstart~=evalind
%         curveind=curveind+1;
%         i;
%         evalstart=evalind;
%         maxlengthind=max(lengthind,maxlengthind);
%         lengthind=0;
%     end
%     datamateq(1,i)=xeqcurve(3,i); % This is the parameter that is varied.
%     datamateq(2,i)=xeqcurve(1,i); % This is the dependent axis of the bifurcation plot.  The one you wish to plot
%     datamateq(3,i)=evalind;
%     datamateq(4,i)=curveind;
%
%     lengthind=lengthind+1;
% end
%
% maxlengthind=max(maxlengthind,lengthind);
%
% curveindeq=curveind;
%
% for i=1:curveindeq
%     index=find(datamateq(4,:)==i);
%     eval(['curve' num2str(i) 'eq' '=datamateq(1:3,index);']);
% end
%
% for i=1:curveindeq
%     stability=eval(['curve' num2str(i) 'eq(3,1)']);
%     if stability==0
%         plotsty='-';
%     else
%         plotsty=':';
%     end
%
%     plotcolor='k';
%
%     plotstr=strcat(plotcolor,plotsty);
%
%     plot(eval(['curve' num2str(i) 'eq(1,:)']),eval(['curve' num2str(i) 'eq(2,:)']),plotstr,'Linewidth',lw)
%     hold on
% end
%
% %add hopf plot
% %start index at 1 to plot G vs beta
% %start index at 2 to plot L vs beta
% ndim = 1;
% M = max(xlc(1:ndim:end-2,:));
% m = min(xlc(1:ndim:end-2,:));
% plot(xlc(end,:),M,xlc(end,:),m,'LineWidth',lw,'Color',highcontrast(4,:))
% axis([0 0.6 0 1.5])
% % xlabel({'$\beta$'},'FontSize',fontsize,'color',[0 0 0],'Interpreter','latex')
% % ylabel({'$G$'},'FontSize',fontsize,'color',[0 0 0],'Interpreter','latex')
% % title('Single GTPase Tension','fontsize',fontsizevar,'color',[0 0 0])
% grid()
% set(gca,'YColor',[0 0 0],'XColor',[0 0 0]);
















%extract HOM-NS point..
%hompts
% jacobian_calc; %calc jacobian symbolically, JFunc(...)
% %@(G,G_T,G_h,L,alpha,b,bta,ell0,epsilon,gama,n,phi) (order of arguments! careful!)
% for j = 1:size(xlc2,2)
%   if trace(JFunc()) == 0
%
%   end

%other homoclinic....
% xtmp = x2(1:2,s2(3).index)
% pvec(aps) = x2(3:4,s2(3).index)
% %initial
% [x0,v0] = init_BT_Hom(syshandle,xtmp,s2(3),pvec,aps,50,4,0.00001,0.0001,[0 1 1]);
% opt=contset;
% opt = contset(opt,'InitStepsize',0.01);
% opt=contset(opt,'MaxStepsize',0.01);
% opt = contset(opt,'MaxNumPoints',3000);
% % opt=contset(opt,'MaxNewtonIters',10);
% [xlc3,vlc3,slc3,hlc3,flc3] = cont(@homoclinic,x0,v0,opt);
% hompts2 = xlc3(end-5:end-4,:);
% plot(hompts2(1,:),hompts2(2,:),'color',purple,'linewidth',1)
%
%
% %this works then crashes
% [x0,v0]=init_GH_LPC(syshandle,x2b,pvec,s2b(3),aps,30,4,0.001);
% opt = contset;
% opt = contset(opt,'MaxNumPoints',500);
% opt=contset(opt,'Singularities',1);
% opt=contset(opt,'MaxStepsize',0.1);
% [x6,v6,s6,h6,f6]=cont(@limitpointcycle,x0,v0,opt); %takes 2185 seconds
% % [x5,v5,s5,h5,f5] = cont(x5,v5,s5,h5,f5,cds); %extend
% plot(x6(end-1,:),x6(end,:),'color',green,'linewidth',1)
% for j = 2:size(x6,2)-1 %don't do the whole curve becuase the starting point and end point shouldn't be flagged
%   if ~isempty(find([s6.index]==j))
%     plot(x6(end-1,j),x6(end,j),'r-',color)
%   end
% end

% make phase plane at SNIC-H pt
% p.b = hompts(2,end);
% p.beta = hompts(1,end);
% p
% GTPasetension;
