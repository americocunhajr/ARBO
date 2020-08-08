
% Plotting Notes
% -----------------------------------------------------------
    % Plots are on a week time scale. Divide 'time' vector by 7. Multiply
        % "day^-1 unit" vectors by 7.
% -----------------------------------------------------------


% Plot number of susceptible humans - SH
% % ...........................................................
% 
% gname  = 'SEIR_SEI_SH';                                                 % file name
% 
% figure('Name',gname,'NumberTitle','off')                                % figure name setting
% 
% fig = plot(time/7,SH,'-b');                                             % plot
%      
%     title('','FontSize',20,'FontName','Helvetica');                     % figure title
%     
%     set(gcf,'color','white');                                           % figure background color
%     
%     xlabel('time (weeks)','FontSize',20,'FontName','Helvetica');        % x axis label and settings
%     ylabel('number of people','FontSize',20,'FontName','Helvetica');    % y axis label and settings
% 
%     xlim([t0/7 t1/7]);                                                  % x axis limit value
%     ylim('auto');                                                       % y axis limit value
%     %axis auto                                                           % axis auto adjust
%     %axis([0.99 52 205.8e6 206.0e6]);                                    % axis manual adjust
%     
%     
%     set(gca,'Box','on');                                                % box around graph 
%     set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);                   % color of the box outline, axis and axis labels
%     set(gca,'TickDir','out','TickLength',[.02 .02]);                    % tick settings
%     set(gca,'XMinorTick','on','YMinorTick','on');                       % Intermediary ticks 
%     set(gca,'FontSize',17,'FontName','Helvetica');                      % tick labels size and font
%     set(gca,'XGrid','off','YGrid','on');                                % grid
%     
%     set(fig,'LineWidth',1.5);                                           % plot line width
%     set(fig,'MarkerSize',2.0);                                          % plot marker size
%     set(fig,'MarkerFaceColor','w');                                     % plot marker fill color
%     set(fig,'MarkerEdgeColor','k');                                     % plot marker edge color
%     
%     % Extra settings
%     %----------------------------------------------------------------------
%         % Settings for y axis tick labels and order of magnitude
%             set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
%             set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes
%             yl = get(gca,'ylim');
%             set(gca,'yTick',linspace(yl(1),yl(2),3))                    % setting number of tick labels to display          
%             BD = 3;                                                     % number of significant figures before the decimal point for the highest tick label
%             OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
%             ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
%             % Formating new tick labels
%                 nyt=cell(size(ryt));
%                 for j=1:length(ryt)
%                     nyt{j}=sprintf('% 1.1f',ryt(j));
%                         % example: '% W.Xf' displays fixed-point notation with X
%                             % digits after the decimal point, minimum of W characters.
%                             % The space after the percent inserts a space before the
%                             % displayed value, giving the same size to + and - numbers. 
%                 end
%             set(gca,'yticklabel',nyt);                                  % setting tick labels
%             % Placing order of magnitude
%                 fs = get(gca,'fontsize');
%                 set(gca,'units','normalized');
%                 xl = xlim;
%                 text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');
%     %----------------------------------------------------------------------
% 
%     
%     saveas(gcf,gname,'epsc2');                                          % saving file
% 
%     %close;                                                             % close figure
% % ...........................................................


% Plot number of incubating humans - EH
% % ...........................................................
% 
% gname  = 'SEIR_SEI_EH';                                                 % file name
% 
% figure('Name',gname,'NumberTitle','off')                                % figure name setting
% 
% fig = plot(time/7,EH,'-b');                                             % plot
%      
%     title('','FontSize',20,'FontName','Helvetica');                     % figure title
%     
%     set(gcf,'color','white');                                           % figure background color
%     
%     xlabel('time (weeks)','FontSize',20,'FontName','Helvetica');        % x axis label and settings
%     ylabel('number of people','FontSize',20,'FontName','Helvetica');    % y axis label and settings
% 
%     xlim([t0/7 t1/7]);                                                  % x axis limit value
%     ylim('auto');                                                       % y axis limit value
%     %axis auto                                                           % axis auto adjust
%     %axis([0.99 52 0 15000]);                                            % axis manual adjust
%     
%     
%     set(gca,'Box','on');                                                % box around graph 
%     set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);                   % color of the box outline
%     set(gca,'TickDir','out','TickLength',[.02 .02]);                    % tick settings
%     set(gca,'XMinorTick','on','YMinorTick','on');                       % Intermediary ticks 
%     set(gca,'FontSize',15,'FontName','Helvetica');                      % tick labels size and font
%     set(gca,'XGrid','off','YGrid','on');                                % grid
%     
%     set(fig,'LineWidth',1.5);                                           % plot line width
%     set(fig,'MarkerSize',2.0);                                          % plot marker size
%     set(fig,'MarkerFaceColor','w');                                     % plot marker fill color
%     set(fig,'MarkerEdgeColor','k');                                     % plot marker edge color
% 
%     % Extra settings
%     %----------------------------------------------------------------------
%         % Settings for y axis tick labels and order of magnitude
%             set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
%             set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes
%             yl = get(gca,'ylim');
%             set(gca,'yTick',linspace(yl(1),yl(2),4))                    % setting number of tick labels to display          
%             BD = 5;                                                     % number of significant figures before the decimal point for the highest tick label
%             OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
%             ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
%             % Formating new tick labels
%                 nyt=cell(size(ryt));
%                 for j=1:length(ryt)
%                     nyt{j}=sprintf('% 1.0f',ryt(j));
%                         % example: '% W.Xf' displays fixed-point notation with X
%                             % digits after the decimal point, minimum of W characters.
%                             % The space after the percent inserts a space before the
%                             % displayed value, giving the same size to + and - numbers. 
%                 end
%             set(gca,'yticklabel',nyt);                                  % setting tick labels
%             % Placing order of magnitude
%                  fs = get(gca,'fontsize');
%                 set(gca,'units','normalized');
%                 xl = xlim;
%                 text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');
%     %----------------------------------------------------------------------
% 
% 
%     saveas(gcf,gname,'epsc2');                                          % saving file
%  
%     %close;                                                             % close figure
% % ........................................................... 


% Plot number of infectious humans - IH
% ...........................................................

gname  = 'SEIR_SEI_IH';                                                 % file name

figure('Name',gname,'NumberTitle','off')                                % figure name setting

fig = plot(time/7,IH,'-b');                                             % plot
     
    title('','FontSize',25,'FontName','Helvetica');                     % figure title
    
    set(gcf,'color','white');                                           % figure background color
    
    xlabel('time (weeks)','FontSize',28,'FontName','Helvetica');        % x axis label and settings
    ylabel('number of people','FontSize',28,'FontName','Helvetica');    % y axis label and settings

    xlim([t0/7 t1/7]);                                                  % x axis limit value
    ylim([0 300e3]);                                                       % y axis limit value
    %axis auto                                                           % axis auto adjust
    %axis([0.99 52 0 15000]);                                            % axis manual adjust
    
    
    set(gca,'Box','on');                                                % box around graph 
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);                   % color of the box outline
    set(gca,'TickDir','out','TickLength',[.02 .02]);                    % tick settings
    set(gca,'XMinorTick','on','YMinorTick','on');                       % Intermediary ticks 
    set(gca,'FontSize',25,'FontName','Helvetica');                      % tick labels size and font
    set(gca,'XGrid','off','YGrid','on');                                % grid
    
    set(fig,'LineWidth',1.5);                                           % plot line width
    set(fig,'MarkerSize',2.0);                                          % plot marker size
    set(fig,'MarkerFaceColor','w');                                     % plot marker fill color
    set(fig,'MarkerEdgeColor','k');                                     % plot marker edge color
    
    
    % Extra settings
    %----------------------------------------------------------------------
        % Settings for y axis tick labels and order of magnitude
            set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
            set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes
            yl = get(gca,'ylim');
            set(gca,'yTick',linspace(yl(1),yl(2),4))                    % setting number of tick labels to display          
            BD = 3;                                                     % number of significant figures before the decimal point for the highest tick label
            OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
            ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
            % Formating new tick labels
                nyt=cell(size(ryt));
                for j=1:length(ryt)
                    nyt{j}=sprintf('% 1.0f',ryt(j));
                        % example: '% W.Xf' displays fixed-point notation with X
                            % digits after the decimal point, minimum of W characters.
                            % The space after the percent inserts a space before the
                            % displayed value, giving the same size to + and - numbers. 
                end
            set(gca,'yticklabel',nyt);                                  % setting tick labels
            % Placing order of magnitude
                fs = get(gca,'fontsize');
                set(gca,'units','normalized');
                xl = xlim;
                text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');
    %----------------------------------------------------------------------
    
    
    %print(gcf,gname,'-dpdf','-r300','-bestfit');
    
    %close;                                                             % close figure
% ...........................................................


% Plot number of recovered humans - RH
% % ...........................................................
% 
% gname  = 'SEIR_SEI_RH';                                                 % file name
% 
% figure('Name',gname,'NumberTitle','off')                                % figure name setting
% 
% fig = plot(time/7,RH,'-b');                                             % plot
%      
%     title('','FontSize',20,'FontName','Helvetica');                     % figure title
%     
%     set(gcf,'color','white');                                           % figure background color
%     
%     xlabel('time (weeks)','FontSize',20,'FontName','Helvetica');        % x axis label and settings
%     ylabel('number of people','FontSize',20,'FontName','Helvetica');    % y axis label and settings
% 
%     xlim([t0/7 t1/7]);                                                  % x axis limit value
%     ylim('auto');                                                       % y axis limit value
%     %axis auto                                                           % axis auto adjust
%     %axis([0.99 52 0 150000]);                                           % axis manual adjust
%     
%     
%     set(gca,'Box','on');                                                % box around graph 
%     set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);                   % color of the box outline
%     set(gca,'TickDir','out','TickLength',[.02 .02]);                    % tick settings
%     set(gca,'XMinorTick','on','YMinorTick','on');                       % Intermediary ticks 
%     set(gca,'FontSize',17,'FontName','Helvetica');                      % tick labels size and font
%     set(gca,'XGrid','off','YGrid','on');                                % grid
%     
%     set(fig,'LineWidth',1.5);                                           % plot line width
%     set(fig,'MarkerSize',2.0);                                          % plot marker size
%     set(fig,'MarkerFaceColor','w');                                     % plot marker fill color
%     set(fig,'MarkerEdgeColor','k');                                     % plot marker edge color
%     
%     
%     % Extra settings
%     %----------------------------------------------------------------------
%         % Settings for y axis tick labels and order of magnitude
%             set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
%             set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes
%             yl = get(gca,'ylim');
%             set(gca,'yTick',linspace(yl(1),yl(2),4))                    % setting number of tick labels to display          
%             BD = 3;                                                     % # of SF before the point in highest tick label (exception: if highest=1 use 0)
%             OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
%             ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
%             % Formating new tick labels
%                 nyt=cell(size(ryt));
%                 for j=1:length(ryt)
%                     nyt{j}=sprintf('% 1.0f',ryt(j));
%                         % example: '% W.Xf' displays fixed-point notation with X
%                             % digits after the decimal point, minimum of W characters.
%                             % The space after the percent inserts a space before the
%                             % displayed value, giving the same size to + and - numbers. 
%                 end
%             set(gca,'yticklabel',nyt);                                  % setting tick labels
%             % Placing order of magnitude
%                 fs = get(gca,'fontsize');
%                 set(gca,'units','normalized');
%                 xl = xlim;
%                 text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');
%     %----------------------------------------------------------------------    
% 
%     
%     print(gcf,gname,'-dpdf','-r300','-bestfit');                                          % saving file
% 
%     %close;                                                             % close figure
% % ...........................................................


% Plot proportion of susceptible vectors - SV  
% % ...........................................................
% 
% gname  = 'SEIR_SEI_SV';                                                 % file name
% 
% figure('Name',gname,'NumberTitle','off')                                % figure name setting
% 
% fig = plot(time/7,SV,'-b');                                             % plot
%      
%     title('','FontSize',20,'FontName','Helvetica');                     % figure title
%     
%     set(gcf,'color','white');                                           % figure background color
%     
%     xlabel('time (weeks)','FontSize',20,'FontName','Helvetica');        % x axis label and settings
%     ylabel('proportion of vectors','FontSize',20,'FontName','Helvetica');    % y axis label and settings
% 
%     xlim([t0/7 t1/7]);                                                  % x axis limit value
%     ylim('auto');                                                       % y axis limit value
%     retriveYlim = get(gca,'ylim');                                      % retrive y axis limit for setting a single limit
%     set(gca,'ylim',[retriveYlim(1) 1]);                                 % setting a single limit
%     %axis auto                                                           % axis auto adjust
%     %axis([0.99 52 0.99993 1]);                                          % axis manual adjust
%     
%     
%     set(gca,'Box','on');                                                % box around graph 
%     set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);                   % color of the box outline
%     set(gca,'TickDir','out','TickLength',[.02 .02]);                    % tick settings
%     set(gca,'XMinorTick','on','YMinorTick','on');                       % Intermediary ticks 
%     set(gca,'FontSize',17,'FontName','Helvetica');                      % tick labels size and font
%     set(gca,'XGrid','off','YGrid','on');                                % grid
%     
%     set(fig,'LineWidth',1.5);                                           % plot line width
%     set(fig,'MarkerSize',2.0);                                          % plot marker size
%     set(fig,'MarkerFaceColor','w');                                     % plot marker fill color
%     set(fig,'MarkerEdgeColor','k');                                     % plot marker edge color
%     
%     % Extra settings
%     %----------------------------------------------------------------------
%         % Settings for y axis tick labels and order of magnitude
%             set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
%             set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes
%             yl = get(gca,'ylim');
%             set(gca,'yTick',linspace(yl(1),yl(2),4))                    % setting number of tick labels to display          
%             BD = 0;                                                     % # of SF before the point in highest tick label (exception: if highest=1 use 0) 
%             OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
%             ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
%             % Formating new tick labels
%                 nyt=cell(size(ryt));
%                 for j=1:length(ryt)
%                     nyt{j}=sprintf('% 1.5f',ryt(j));
%                         % example: '% W.Xf' displays fixed-point notation with X
%                             % digits after the decimal point, minimum of W characters.
%                             % The space after the percent inserts a space before the
%                             % displayed value, giving the same size to + and - numbers. 
%                 end
%             set(gca,'yticklabel',nyt);                                  % setting tick labels
%             % Placing order of magnitude
%                 fs = get(gca,'fontsize');
%                 set(gca,'units','normalized');
%                 xl = xlim;
%                 text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');        
%     %----------------------------------------------------------------------
% 
%     
%     print(gcf,gname,'-dpdf','-r300','-bestfit');                                          % saving file
% 
%     %close;                                                             % close figure
% % ...........................................................


% Plot proportion of incubating vectors - EV  
% % ...........................................................
% 
% gname  = 'SEIR_SEI_EV';                                                 % file name
% 
% figure('Name',gname,'NumberTitle','off')                                % figure name setting
% 
% fig = plot(time/7,EV,'-b');                                             % plot
%      
%     title('','FontSize',20,'FontName','Helvetica');                     % figure title
%     
%     set(gcf,'color','white');                                           % figure background color
%     
%     xlabel('time (weeks)','FontSize',20,'FontName','Helvetica');        % x axis label and settings
%     ylabel('proportion of vectors','FontSize',20,'FontName','Helvetica');    % y axis label and settings
% 
%     xlim([t0/7 t1/7]);                                                  % x axis limit value
%     ylim('auto');                                                       % y axis limit value
%     %axis auto                                                           % axis auto adjust
%     %axis([0.99 52 0 0.00015]);                                          % axis manual adjust
%     
%     
%     set(gca,'Box','on');                                                % box around graph 
%     set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);                   % color of the box outline
%     set(gca,'TickDir','out','TickLength',[.02 .02]);                    % tick settings
%     set(gca,'XMinorTick','on','YMinorTick','on');                       % Intermediary ticks 
%     set(gca,'FontSize',17,'FontName','Helvetica');                      % tick labels size and font
%     set(gca,'XGrid','off','YGrid','on');                                % grid
%     
%     set(fig,'LineWidth',1.5);                                           % plot line width
%     set(fig,'MarkerSize',2.0);                                          % plot marker size
%     set(fig,'MarkerFaceColor','w');                                     % plot marker fill color
%     set(fig,'MarkerEdgeColor','k');                                     % plot marker edge color
%     
%     % Extra settings
%     %----------------------------------------------------------------------
%         % Settings for y axis tick labels and order of magnitude
%             set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
%             set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes
%             yl = get(gca,'ylim');
%             set(gca,'yTick',linspace(yl(1),yl(2),4))                    % setting number of tick labels to display          
%             BD = 1;                                                     % # of SF before the point in highest tick label (exception: if highest=1 use 0)
%             OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
%             ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
%             % Formating new tick labels
%                 nyt=cell(size(ryt));
%                 for j=1:length(ryt)
%                     nyt{j}=sprintf('% 1.1f',ryt(j));
%                         % example: '% W.Xf' displays fixed-point notation with X
%                             % digits after the decimal point, minimum of W characters.
%                             % The space after the percent inserts a space before the
%                             % displayed value, giving the same size to + and - numbers. 
%                 end
%             set(gca,'yticklabel',nyt);                                  % setting tick labels
%             % Placing order of magnitude
%                 fs = get(gca,'fontsize');
%                 set(gca,'units','normalized');
%                 xl = xlim;
%                 text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');
%     %----------------------------------------------------------------------
% 
%     
%     print(gcf,gname,'-dpdf','-r300','-bestfit');                                          % saving file
% 
%     %close;                                                             % close figure
% % ...........................................................


% Plot proportion of infectious vectors - IV
% % ...........................................................
% 
% gname  = 'SEIR_SEI_IV';                                                 % file name   
% 
% figure('Name',gname,'NumberTitle','off')                                % figure name setting
% 
% fig = plot(time/7,IV,'-b');                                             % plot
%      
%     title('','FontSize',20,'FontName','Helvetica');                     % figure title
%     
%     set(gcf,'color','white');                                           % figure background color
%     
%     xlabel('time (weeks)','FontSize',20,'FontName','Helvetica');        % x axis label and settings
%     ylabel('proportion of vectors','FontSize',20,'FontName','Helvetica');    % y axis label and settings
% 
%     xlim([t0/7 t1/7]);                                                  % x axis limit value
%     ylim('auto');                                                       % y axis limit value
%     %axis auto                                                           % axis auto adjust
%     %axis([0.99 52 0 0.00015]);                                          % axis manual adjust
%     
%     
%     set(gca,'Box','on');                                                % box around graph 
%     set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);                   % color of the box outline
%     set(gca,'TickDir','out','TickLength',[.02 .02]);                    % tick settings
%     set(gca,'XMinorTick','on','YMinorTick','on');                       % Intermediary ticks 
%     set(gca,'FontSize',17,'FontName','Helvetica');                      % tick labels size and font
%     set(gca,'XGrid','off','YGrid','on');                                % grid
%     
%     set(fig,'LineWidth',1.5);                                           % plot line width
%     set(fig,'MarkerSize',2.0);                                          % plot marker size
%     set(fig,'MarkerFaceColor','w');                                     % plot marker fill color
%     set(fig,'MarkerEdgeColor','k');                                     % plot marker edge color
%     
%     % Extra settings
%     %----------------------------------------------------------------------
%         % Settings for y axis tick labels and order of magnitude
%             set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
%             set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes        
%             yl = get(gca,'ylim');
%             set(gca,'yTick',linspace(yl(1),yl(2),4))                    % setting number of tick labels to display          
%             BD = 1;                                                     % # of SF before the point in highest tick label (exception: if highest=1 use 0)
%             OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
%             ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
%             % Formating new tick labels
%                 nyt=cell(size(ryt));
%                 for j=1:length(ryt)
%                     nyt{j}=sprintf('% 1.1f',ryt(j));
%                         % example: '% W.Xf' displays fixed-point notation with X
%                             % digits after the decimal point, minimum of W characters.
%                             % The space after the percent inserts a space before the
%                             % displayed value, giving the same size to + and - numbers. 
%                 end
%             set(gca,'yticklabel',nyt);                                  % setting tick labels
%             % Placing order of magnitude
%                 fs = get(gca,'fontsize');
%                 set(gca,'units','normalized');
%                 xl = xlim;
%                 text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');
%     %----------------------------------------------------------------------
%     
%     
%     print(gcf,gname,'-dpdf','-r300','-bestfit');                                          % saving file
% 
%     %close;                                                             % close figure
% % ...........................................................


% Plot cumulative number of infectious humans - C
% ...........................................................

gname  = 'SEIR_SEI_C';                                                  % file name

figure('Name',gname,'NumberTitle','off')                                % figure name setting


fig = plot(time/7,C,'-b');                                              % plot
     
    title('','FontSize',25,'FontName','Helvetica');                     % figure title
    
    set(gcf,'color','white');                                           % figure background color
    
    xlabel('time (weeks)','FontSize',28,'FontName','Helvetica');        % x axis label and settings
    ylabel('number of people','FontSize',28,'FontName','Helvetica');    % y axis label and settings

    xlim([t0/7 t1/7]);                                                  % x axis limit value
    ylim([0 300e3]);                                                       % y axis limit value
    %axis auto                                                           % axis auto adjust
    %axis([0.99 52 0 150000]);                                           % axis manual adjust
    
    
    set(gca,'Box','on');                                                % box around graph 
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);                   % color of the box outline
    set(gca,'TickDir','out','TickLength',[.02 .02]);                    % tick settings
    set(gca,'XMinorTick','on','YMinorTick','on');                       % Intermediary ticks 
    set(gca,'FontSize',25,'FontName','Helvetica');                      % tick labels size and font
    set(gca,'XGrid','off','YGrid','on');                                % grid
    
    set(fig,'LineWidth',1.5);                                           % plot line width
    set(fig,'MarkerSize',2.0);                                          % plot marker size
    set(fig,'MarkerFaceColor','w');                                     % plot marker fill color
    set(fig,'MarkerEdgeColor','k');                                     % plot marker edge color
    
    % Extra settings
    %----------------------------------------------------------------------
        % Scatter plot
          hold on
          scatter(day/7,CData,'LineWidth',1.2,'MarkerEdgeColor','red')  % data plot
          hold off
          
        % Settings for y axis tick labels and order of magnitude
            set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
            set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes
            yl = get(gca,'ylim');                                       
            set(gca,'yTick',linspace(yl(1),yl(2),4))                    % setting number of tick labels to display          
            BD = 3;                                                     % # of SF before the point in highest tick label (exception: if highest=1 use 0)
            OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
            ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
            % Formating new tick labels
                nyt=cell(size(ryt));
                for j=1:length(ryt)
                    nyt{j}=sprintf('% 1.0f',ryt(j));
                        % example: '% W.Xf' displays fixed-point notation with X
                            % digits after the decimal point, minimum of W characters.
                            % The space after the percent inserts a space before the
                            % displayed value, giving the same size to + and - numbers. 
                end
            set(gca,'yticklabel',nyt);                                  % setting tick labels
            % Placing order of magnitude
                fs = get(gca,'fontsize');
                set(gca,'units','normalized');
                xl = xlim;
                text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');
            
            set(gca,'xTick',[0 10 20 30 40 50]);
    %----------------------------------------------------------------------
    
    
     print(gcf,gname,'-dpdf','-r300','-bestfit');                       % saving file
     
    %close;                                                             % close figure
% ...........................................................


% Stem plot new number of cases - NewCases
% ...........................................................

gname  = 'SEIR_SEI_NewCases';                                           % file name

figure('Name',gname,'NumberTitle','off')                                % figure name setting


fig1 = stem(NewCasesData(1:end-2));                                              % data plot
    
hold on
fig2 = stem(NewCases(1:end));                                                  % system response
hold off


    title('','FontSize',25,'FontName','Helvetica');                     % figure title

    leg1 = 'data';                                                      % plot 1 legend     
    leg2 = 'model';                                                     % plot 2 legend
    leg = legend(leg1,leg2);                        
    set(leg,'FontName','Helvetica','FontSize',25);                      % legend font and size

    set(gcf,'color','white');                                           % figure background color
    
    xlabel('time (weeks)','FontSize',28,'FontName','Helvetica');    % x axis label and settings
    ylabel('number of people','FontSize',28,'FontName','Helvetica');               % y axis label and settings

    xlim([0.99 t1/7]);                                                  % x axis limit value
    ylim([0 25000]);                                                       % y axis limit value
    %axis auto                                                           % axis auto adjust
    %axis([0.99 52 0 150000]);                                           % axis manual adjust
    
    
    set(gca,'Box','on');                                                % box around graph 
    set(gca,'XColor',[.3 .3 .3],'YColor',[.3 .3 .3]);                   % color of the box outline
    set(gca,'TickDir','out','TickLength',[.02 .02]);                    % tick settings
    set(gca,'XMinorTick','on','YMinorTick','on');                       % Intermediary ticks 
    set(gca,'FontSize',25,'FontName','Helvetica');                      % tick labels size and font
    set(gca,'XGrid','off','YGrid','on');                                % grid
    
    set(fig1,'LineStyle','--','LineWidth',1.4);                          % plot1 line style and width
    set(fig2,'LineStyle','--','LineWidth',1.4);                          % plot2 line style and width
    set(fig1,'Marker','o','MarkerSize',8.0);                            % plot1 marker symbol and size
    set(fig2,'Marker','s','MarkerSize',5.0);                            % plot2 marker symbol and size
    set(fig1,'Color','r');                                              % plot1 line color
    set(fig2,'Color','b');                                              % plot2 line color
    set(fig1,'MarkerFaceColor','none');                                    % plot1 marker fill color
    set(fig2,'MarkerFaceColor','b');                                    % plot2 marker fill color
    set(fig1,'MarkerEdgeColor','r');                                    % plot1 marker edge color
    set(fig2,'MarkerEdgeColor','b');                                    % plot2 marker edge color
    
    % Extra settings
    %----------------------------------------------------------------------
        set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
        set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes
        
        % Settings for y axis tick labels and order of magnitude
            set(gca,'XTickMode','manual','YTickMode','manual')          % preserve tick values for all figure sizes
            set(gca,'XLimMode','manual','YLimMode','manual')            % preserve axis limits for all figure sizes
            yl = get(gca,'ylim');                                       
            set(gca,'yTick',linspace(yl(1),yl(2),5))                    % setting number of tick labels to display          
            BD = 2;                                                     % # of SF before the point in highest tick label (exception: if highest=1 use 0)
            OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
            ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
            % Formating new tick labels
                nyt=cell(size(ryt));
                for j=1:length(ryt)
                    nyt{j}=sprintf('% 1.0f',ryt(j));
                        % example: '% W.Xf' displays fixed-point notation with X
                            % digits after the decimal point, minimum of W characters.
                            % The space after the percent inserts a space before the
                            % displayed value, giving the same size to + and - numbers. 
                end
            set(gca,'yticklabel',nyt);                                  % setting tick labels
            % Placing order of magnitude
                fs = get(gca,'fontsize');
                set(gca,'units','normalized');
                xl = xlim;
                text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');
                
%         % Displaying curve information
%         Peak   = max(NewCases(1:end-2));                                                          % maximum value of NewCases
%         Diff   = max(NewCasesData(1:end-2)) - Peak;                                               % difference of peak and maximum value of NewCasesData
%         tPeak  = time( ( 7*find(NewCases(1:end-2)==max(NewCases(1:end-2))) - time(1) )/dt +1 )/7; % week of the peak on IH
%         tDiff  = find(NewCasesData(1:end-2)==max(NewCasesData(1:end-2))) - tPeak;        % week difference of peak and maximum value of data confirmed_IH
% 
%         text(18,22500,['Peak = ',num2str(Peak)],'Fontsize',10);
%         text(18,21000,['Diff = ',num2str(Diff)],'Fontsize',10);
%         text(18,19500,['t_{Peak} = ',num2str(tPeak)],'Fontsize',10);
%         text(18,18000,['t_{Diff} = ',num2str(tDiff)],'Fontsize',10);
% %bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 IV0 (case4)
%         text(28,15500,['N = ',num2str(N)],'Fontsize',10);
%         text(28,14000,['TaH = ',num2str(TaH)],'Fontsize',10);
%         text(28,12500,['TaV = ',num2str(TaV)],'Fontsize',10);
%         text(28,11000,['TyH = ',num2str(TyH)],'Fontsize',10);
%         text(28,9500,['TdV = ',num2str(TdV)],'Fontsize',10);
%         text(28,8000,['bH = ',num2str(bH)],'Fontsize',10);
%         text(28,6500,['bV = ',num2str(bV)],'Fontsize',10);
% 
%         set(gca,'Box','off');
%         text(42,15500,['SH0 = ',num2str(SH0)],'Fontsize',10);
%         text(42,14000,['EH0 = ',num2str(EH0)],'Fontsize',10);
%         text(42,12500,['IH0 = ',num2str(IH0)],'Fontsize',10);
%         text(42,11000,['RH0 = ',num2str(RH0)],'Fontsize',10);
%         text(42,9500,['C0 = ',num2str(C0)],'Fontsize',10);
%         text(42,8000,['SV0 = ',num2str(SV0)],'Fontsize',10);
%         text(42,6500,['IV0 = ',num2str(IV0)],'Fontsize',10);
%         text(42,5000,['EV0 = ',num2str(EV0)],'Fontsize',10);             

    %----------------------------------------------------------------------
    
    print(gcf,gname,'-dpdf','-r300','-bestfit');
     
    %close
% ...........................................................     
    