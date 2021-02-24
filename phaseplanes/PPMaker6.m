function PPmaker5
%UPDATE BY CZ JAN 28 2021, loop through xppaut files, make phaseplanes
%UPDATE BY CZ FEB 1 2021, to match colors in the rest of the paper
%


% This function opens a file browser for the current pwd listing all .dat
% files and plots the solution in the chosen file using the file name to
% determine the solution type.

% district number - not used for now - I assume the user is only going to
% be adding solutions from a single parameter set at a time. The current
% folder name should match the prefixes in the folder. That is, folder name
% District6 should contain .dat files and a parameter .m named
% District6-parameters.m and District6-blah-blah.dat.
  close all
  forExport=1;

  districts = {'d01','d02','d03','d04','d05','d06','d07','d08','d09','d10','d11','d12'};
  titles = {'i','ii','iii','iv','v','vi','vii','viii','ix','x','xi','xii'}
  initialPath=pwd;       % get name of current path

  for i=1:length(districts)
    strPWD=districts{i};    % let user pick a folder/district containing the .dat files
    cd(strPWD)          % change to that folder/district
    % extract the folder name from the full path
    pathParts=regexp(strPWD,filesep,'split');
    folder=pathParts{end};
    % import parameters.dat which contained the parameter values for the
    % current folder/district, beta first, b second
    params=importdata('parameters.dat');
    beta=params(1);
    b=params(2);

    dirInfo=dir('.');           % get a list of files in the folder
    numFiles=size(dirInfo,1);   % get count number of files in folder

    bettercolors;

    figure(1)
    clf
    % title('Do you want to export to eps when done? y/n')
    % k=waitforbuttonpress;if get(gcf,'CurrentCharacter')=='y', forExport=1; end

    % strTitle=sprintf('Phase plane for beta=%0.5f, b=%0.5f.',beta,b)
    % title(strTitle)

    % strTitle=sprintf('Phase plane for \\beta = %0.5f, b = %0.5f',beta,b)
    % title(strTitle,'fontweight','normal')
    %get rid of title since parameters are in the filename
    % xlabel('$G$','Interpreter','latex')
    % ylabel('$L$','Interpreter','latex')
    set(gca,'YColor','k')
    set(gca,'XColor','k')
    set(gca,'box','on')

    hold on

    % loop over all files
    plotFromDatFiles(numFiles,dirInfo,folder)

    % plot some other solutions for context (interactive)

    p.alpha=10;
    p.gamma=1.5;
    p.n=4;
    p.G_T=2;
    p.phi=0.75;
    p.G_h=0.3;
    p.epsilon=0.1;
    p.ell0=1;

    p.beta=beta;
    p.b=b;

    MoreSolutions=0;

    while MoreSolutions

        title('Click to add a solution.')
        x=ginput(1);

        y0 = [x(1),x(2)]
        Time1 = [0,50];
        Time2 = [0,-50];
        [t1,y1] = ode45(@(t,y) odes(t,y,p),Time1,y0);
        [t2,y2] = ode45(@(t,y) odes(t,y,p),Time2,y0);
        y=[flip(y2);y1];
        t=[flip(t2);t1];
        h=plot(y(:,1),y(:,2),'color',grey,'LineWidth',0.5);

        title('Keep this solution? y/n')
        k=waitforbuttonpress;
        if get(gcf,'CurrentCharacter')=='n'
            delete(h)
        else
            k=1;
            str=[folder '-orbit-n%d.dat'];
            strNewDatFileName=sprintf(str,k)
            while isfile(strNewDatFileName)
                k=k+1;
                strNewDatFileName=sprintf(str,k);
            end
            fileID = fopen(strNewDatFileName,'w');
            fprintf(fileID,'%0.5f %0.5f %0.5f\n',[t y(:,1) y(:,2)]');
            fclose(fileID);

        end

        title('Add more solutions? y/n')
        k=waitforbuttonpress;if get(gcf,'CurrentCharacter')=='n', MoreSolutions=0; end

    end

    plotFromDatFiles(numFiles,dirInfo,folder)

    x0=5;y0=5;

    cd(initialPath)

    if forExport
        happyZoom=1;

        while ~happyZoom
            title('Click on two diagonal corners to set the zoom.')
            coords=ginput(2);
            xmin=min(coords(:,1));
            ymin=min(coords(:,2));
            xmax=max(coords(:,1));
            ymax=max(coords(:,2));
            axis([xmin xmax ymin ymax])
            title('Happy with this zoom? y/n')
            k=waitforbuttonpress;if get(gcf,'CurrentCharacter')=='y', happyZoom=1; end
        end

        title('Use built-in zoom, then hit d (done).')
        while ~happyZoom
            k=waitforbuttonpress;
            if get(gcf,'CurrentCharacter')=='d', happyZoom=1; end
        end

        %set final zoom
        xlim([0.25,1])
        xticks(linspace(0.25,1,5))
        ylim([0.25,0.52])
        yticks(linspace(0.25,0.52,5))
        xticklabels([])
        yticklabels([])
        axis square

        title('')
        text(0.9,0.48,titles{i},'HorizontalAlignment','center')
        set(gcf,'Units','inches','Position',[x0 y0 2 2],'PaperPositionMode','auto')
        set(gca,'fontsize',10)
        fileTitle = sprintf('beta_%.5f_b%.5f',beta,b)
        print(1,[folder,'_',fileTitle,'.jpg'],'-djpeg','-r600') %dpi = 600 %format = png
        print(1,[folder,'_',fileTitle,'.eps'],'-depsc','-painters') %color eps

    %     filePathTitle=[folder,'_',fileTitle,'.eps'];
    %     exportgraphics(gcf,filePathTitle)
    %     export_fig(gcf,filePathTitle)
    end
  end
end

function plotFromDatFiles(numFiles,dirInfo,folder)

  bettercolors;

    for k=3:numFiles
        strFilename=dirInfo(k).name;    % extract file name from list
        if strFilename(1:3)==folder     % use only files that start with the same string as the folder itself
            type=strFilename(5:9);      % extract the solution type from the file name (orbit, manif, perio)
            stability=strFilename(11);  % extract stability...
            linewidth=[];
            switch type                 % switch conditions for setting line color, thickness, symbols for steady states
                case 'manif'
                    if stability=='s'
    %                    clr=[0 0.7 0];
                        % clr=[0 1 1];    % cyan
                        clr = [0.2,0.2,0.2];
                        linewidth=0.5;
                        style='-';
                        vis=1;
                    elseif stability=='u'
    %                    clr=[0 0.7 0];
                        % clr=[1 1 0];    % yellow
                        clr = [0.8,0.8,0.8];
                        linewidth=0.5;
                        style='-';
                        vis=1;
                    end
                case 'orbit'
                    clr = grey;
                    linewidth=0.5;
                    style='-';
                    vis=0;
                case 'perio'
                    if stability=='s'
                        % clr=[0 1 0];    % green
                        clr = 1-(1-green)*1.1;
                        linewidth=2;
                        style='-';
                        vis=1;
                    elseif stability=='u'
                        % clr=[0 0 1];    % blue
                        clr = 1-(1-green)*0.8;
                        linewidth=1;
                        style='-';
                        vis=1;
                    end
                case 'stead'
                  if stability=='s'
                      % clr=[1 0 0];    % red
                      clr = red.*0.6;
                      linewidth=2;
                      style='o';
                      vis=1;
                  elseif stability=='u'
    %                  clr=[0 0 1];
                      clr=red;    % black
                      linewidth=1;
                      style='o';
                      vis=1;
                  elseif stability=='n'
                          clr=[0 0 0];
                          linewidth=1;
                          style='d';
                          vis=1;
                  end

            end
            if ~isempty(linewidth)
                Solution=importdata(dirInfo(k).name);
                plot(Solution(:,2),Solution(:,3),style,'color',clr,'linewidth',linewidth,'Markersize',4,'visible',vis)
            end
        end
    end

end

function dy = odes(t,y,p)
  dy = zeros(size(y));
  G = y(1);
  L = y(2);
  dy(1) = dGdt(G,L,p);
  dy(2) = dLdt(G,L,p);
end

function dy = dLdt(G,L,p)
  dy = -p.epsilon*(L-(p.ell0-p.phi*G.^p.n./(p.G_h.^p.n+G.^p.n)));
end

function dy = dGdt(G,L,p)
    dy = (p.b+(p.beta./(1+exp(-p.alpha*(L-(p.ell0-p.phi*G.^p.n./(p.G_h.^p.n+G.^p.n))))))+p.gamma*(G.^p.n)./(1+G.^p.n)).*(p.G_T - G)-G;
end
