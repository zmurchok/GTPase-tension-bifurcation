clear all
close all

bettercolors;

load('Fig1.mat');
load('2par-HSN-to-upper-NCH.mat')

idxrange = size(xlc2,2) + size(x,2);

for idx=1:1:idxrange
  if idx <= size(xlc2,2)
    figure()
    set(gcf,'Units','inches','Position',[5,5,3.25*2,3.25*3/4])
    subplot(1,2,1)

    hold on
    %SN
    sn = plot(x3(end-1,:),x3(end,:),'color',red*0.6,'linewidth',1)
    plot(x3b(end-1,:),x3b(end,:),'color',red*0.6,'linewidth',1)
    plot(x4(end-1,:),x4(end,:),'color',red*0.6,'linewidth',1)
    snich = [0.16509656,0.14084249]
    snic = plot([x4(end-1,1:515),snich(1)],[x4(end,1:515),snich(2)],'color',red*0.5+yellow*0.5,'linewidth',1)
    plot(x4b(end-1,:),x4b(end,:),'color',red*0.6,'linewidth',1)

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
    plot(homns(1),homns(2),'k.')
    plot(snich(1),snich(2),'k.')

    set(gca,'YColor','k','Xcolor','k','box','on','fontsize',10)
    ylabel('$b$','Interpreter','latex')
    xlabel('$\beta$','Interpreter','latex')

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

    scatter(hompts(1,idx),hompts(2,idx),'r')
    axis([0.1    0.25        0.12    0.21])


    subplot(1,2,2)
    plot(xlc2(1:2:end-6,idx),xlc2(2:2:end-6,idx),'color',grey)
    axis([0.2    1.2        0.2    0.6])
    xticks([0.2,0.4,0.6,0.8,1,1.2])
    set(gca,'YColor','k','Xcolor','k','box','on','fontsize',10)
    ylabel('$L$','Interpreter','latex')
    xlabel('$G$','Interpreter','latex')
    print(1,sprintf('imgs/%04d.png',idx),'-dpng','-r300')
    close all
  else
    tempidx = size(x,2) - (idx - size(xlc2,2))
    figure
    set(gcf,'Units','inches','Position',[5,5,3.25*2,3.25*3/4])
    subplot(1,2,1)

    hold on
    %SN
    sn = plot(x3(end-1,:),x3(end,:),'color',red*0.6,'linewidth',1)
    plot(x3b(end-1,:),x3b(end,:),'color',red*0.6,'linewidth',1)
    plot(x4(end-1,:),x4(end,:),'color',red*0.6,'linewidth',1)
    snich = [0.16509656,0.14084249]
    snic = plot([x4(end-1,1:515),snich(1)],[x4(end,1:515),snich(2)],'color',red*0.5+yellow*0.5,'linewidth',1)
    plot(x4b(end-1,:),x4b(end,:),'color',red*0.6,'linewidth',1)

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
    plot(homns(1),homns(2),'k.')
    plot(snich(1),snich(2),'k.')

    set(gca,'YColor','k','Xcolor','k','box','on','fontsize',10)
    ylabel('$b$','Interpreter','latex')
    xlabel('$\beta$','Interpreter','latex')

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


    scatter(x(end-4,tempidx),x(end-3,tempidx),'r')
    axis([0.1    0.25        0.12    0.21])


    subplot(1,2,2)
    plot(x(1:2:end-5,tempidx),x(2:2:end-5,tempidx),'color',grey)
    axis([0.2    1.2        0.2    0.6])
    xticks([0.2,0.4,0.6,0.8,1,1.2])
    set(gca,'YColor','k','Xcolor','k','box','on','fontsize',10)
    ylabel('$L$','Interpreter','latex')
    xlabel('$G$','Interpreter','latex')
    print(1,sprintf('imgs/%04d.png',idx),'-dpng','-r300')
    close all
  end
end
