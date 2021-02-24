%https://personal.sron.nl/~pault/#sec:colour_blindness
bright = [68,119,170; 102,204,238; 34,136,51; 204,187,68; 238,102,119; 170,51,119; 187,187,187]/255;
highcontrast = [0,0,0; 0,68,136; 187,85,102; 221,170,51]/255;
muted = [51,24,136;136,204,238;68,170,153;17,119,51;153,153,51;221,204,119;204,102,119;136,34,85;170,68,153]/255;
iridescent = [254,251,233;252,247,213;245,243,193;234,240,181;221,236,191;208,231,202;194,227,210;181,221,216;168,216,220;155,210,225;141,203,228;129,196,231;123,188,231;126,178,228;136,165,221;147,152,210;155,138,196;157,125,178;154,112,158;144,99,136;128,87,112;104,73,87;70,53,58]/255;

blue = bright(1,:);
cyan = bright(2,:);
green = bright(3,:);
yellow = bright(4,:);
red = bright(5,:);
purple = bright(6,:);
grey = bright(7,:);

set(groot,'defaultAxesColorOrder',bright)

% x = linspace(0,15,200);
%
% width=4;
% height=3;
% x0 = 5;
% y0 = 5;
% fontsize = 12;
% figure('Units','inches','Position',[x0 y0 width height],'PaperPositionMode','auto');
% hold on
% for i=1:10
%   plot(x,sin(x+i/2)*(10-i),'LineWidth',2)%'Color',iridescent(i,:),'LineWidth',2);
% end
%
% xlabel({'$x$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
% ylabel({'$y$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
% set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',fontsize,'FontName','Times')
% set(gca,'Box','on')
% set(gca,'LineWidth',1.5)
% grid
%
% figure('Units','inches','Position',[x0 y0 width height],'PaperPositionMode','auto');
% [X,Y] = meshgrid(-8:.5:8);
% R = sqrt(X.^2 + Y.^2) + eps;
% Z = sin(R)./R;
% surf(X,Y,Z)
% xlabel({'$x$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
% ylabel({'$y$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
% zlabel({'$z$'},'FontUnits','points','Interpreter','latex','FontWeight','normal','FontSize',fontsize,'FontName','Times')
% set(gca,'Units','normalized','FontUnits','points','FontWeight','normal','FontSize',fontsize,'FontName','Times')
% set(gca,'Box','on')
% set(gca,'LineWidth',1.5)
% c = colorbar;
% c.LineWidth=1.5;
% colormap(iridescent)
