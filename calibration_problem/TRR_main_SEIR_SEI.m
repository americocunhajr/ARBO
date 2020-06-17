%% Main file for the Trust-Region-Reflective (IC) lsqcurvefit for the SEIR-SEI ZIKAV Brazil outbreak 
% ----------------------------------------------------------------- 
%  Calculations are made on a day time scale and plotting on a week time
%  scale
%   
%   Procedure:
%   The function lsqcurvefit(fun,inital_guess_param,xdata,ydata) provides
%   values for the varying parameters that best approximate F(p,xdata) to
%   the ydata points on a least-squared-differences sense
%   - fun is the function that supply F(p,xdata)
%   - p is the vector containing the varying parameters.
%   - xdata is the data that correspond to the domain values
%   - F(p,xdata) is the function values that must approximate to ydata
%   - ydata is the fel data
%   - initial_guess_param is the initial guess of p that must be supplied
%   to the algorithm
%
%   When using a ode45 function, a second function must be made that provides
%   F(p,xdata). F(p,xdata) will not be the full result of ode45
%   integration. This second function must evaluate the 'y' vector response 
%   of ode45 and pick the F(p,xdata) values on the 'y' vector that 
%   correspond to the ydata values for comparison.
%   A third function is required to integrate through ode45, the
%   system defining function of the integration.
%   The associated fixed and varying parameters must be specified in each
%   function.
%
% ----------------------------------------------------------------- 
%   ZikaV information
%   
%   N = human population size (number of ppl)
%
%   SH = susceptible humans (number of ppl at t)
%   EH = humans incubating (number of ppl at t)
%   IH = infectious humans (number of ppl at t)
%   RH = recovered humans (number of ppl at t)
%
%   bH = vector-to-human transmission rate (days^-1)
%   aH = incubation human frequency (days^-1)
%   yH = infectious human frequency (days^-1)
%
%   C  = cumulative number of infectious humans (number of ppl at t)
%
%
%   SV = proportion of susceptible vector (adimensional)
%   EV = proportion of vectors incubating (adimensional)
%   IV = proportion of infectious vectors (adimensional)
%
%   bV = human-to-vector transmission rate (days^-1)
%   aV = incubation vector frequency (days^-1)
%   dV = vector lifespan frequency (days^-1)
% 
% ----------------------------------------------------------------- 
%  programmers: Eber Dantas
%               eber.paiva@uerj.br
%
%               Michel Tosin
%               michel.tosin@uerj.br
%
%               Americo Cunha
%               americo.cunha@uerj.br
%
%  last update: July 7, 2018
% -----------------------------------------------------------------


%% Header
% -----------------------------------------------------------
clc
clear all
close all


%% Data - Brazil
% -----------------------------------------------------------------
% Load 2016 data
ProbData = load('2016zika_Prob250117.dat');

ConfData = load('2016zika_Conf250117.dat');

% Data time vector
day = [ProbData(1,:)*7];% 400 450];                        % last values are arbitrary, merely for comparison in the algorithm
%day = ConfData(1,:)*7;

% New cases data
NewCasesData = [ProbData(2,:)];% 206e6-29639 1];           % last values are N-RH0 and total mosquito population (in proportion)
%NewCasesData = ConfData(2,:);

% Cumulative cases data
CData = cumsum(NewCasesData);

% Load 2015 data - Cumulative number of cases
%C2015 = 27127;     % until EW 45
C2015 = 29639;      % until EW 49

%C2015 = 443502;     % inferior PAHO estimation     
%C2015 = 1301140;    % superior PAHO estimation
% -----------------------------------------------------------


%% Physical parameters and initial conditions
% -----------------------------------------------------------      
% Human population size (number of ppl)
    N = 206e6;          

% Intrinsic incubation period (days)
    TaH = 5.9;
        
% Extrinsic incubation period (days)
    TaV = 9.1;
        
% Human infectious period (days)
    TyH = 7.9;
        
% Mosquito lifespan (days)
    TdV = 11.0;        
    
% Rates and inverses
    aH  = 1./TaH;                   % intrinsic incubation ratio (days^-1)
    aV  = 1./TaV;                   % extrinsic incubation ratio (days^-1)
    yH  = 1./TyH;                   % human infectious ratio (days^-1)
    dV  = 1./TdV;                   % mosquito death rate (lifespan inverse) (days^-1)

% Transmission rates
    %bH  = 1/11.3;                   % vector-to-human transmission rate (days^-1)
    bH  = 0.1;
    %bV  = 1/8.6;                    % human-to-vector transmission rate (days^-1)
    bV  = 0.1;

% physical parameters vector
    phys_param = [N bH aH yH bV aV dV];
        

% Initial conditions
    IH0 = 300000;          % initial number of infectious humans
        
    EH0 = IH0;                      % initial number of incubating humans
        
    SH0 = N - IH0 - EH0 - C2015;    % initial number of susceptible humans
        
    RH0 = C2015;                        % initial number of recovered humans
        
    C0  = NewCasesData(1);                      % initial cumulative number of infectious humans

        
    IV0 = 0;                  % initial proportion of infectious vectors
        
    EV0 = IV0;                      % initial proportion of incubating vectors
        
    SV0 = 1-EV0-IV0;            % initial proportion of susceptible vectors
     
    % Initial conditions vector
    IC = [SH0 EH0 IH0 SV0 EV0 IV0 C0];      
        % RH0 is omitted to optimize the procedure 
% -----------------------------------------------------------


%% TRR case configuration
% -----------------------------------------------------------
% TRRparam vector includes the system parameters and initial conditions
% Configuration must include the fixed parameters/IC and the initial guess for the variable TRRparameters. 
% Each case configuration modify the other 2 associated files. 
% Maintain default TRRparameter vector order!

% Default TRRparam vector: [N bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 IV0 C0]
    % RH0 is omitted to optimize the procedure
    
% Choose case
case_name = 'case_config5'; 


% Case configuration specification


if strcmp(case_name,'case_config1')

    fixed_TRRparam = [N bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 C0];  
    inital_guess_TRRparam = [IV0]; % variable_TRRparam
    
    
elseif strcmp(case_name,'case_config2')

     fixed_TRRparam = [N bH aH yH bV aV dV SH0 EH0 IH0 SV0 C0];  
     inital_guess_TRRparam = [EV0 IV0]; % variable_TRRparam
 
elseif strcmp(case_name,'case_config3')
    
    fixed_TRRparam = [N aH yH aV dV SH0 EH0 IH0 SV0 C0];  
    inital_guess_TRRparam = [bH bV EV0 IV0]; % variable_TRRparam

elseif strcmp(case_name,'case_config4')
    
    fixed_TRRparam = [N C0];  
    inital_guess_TRRparam = [bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 IV0]; % variable_TRRparam

elseif strcmp(case_name,'case_config5')
    
    % Case with Sums.
    fixed_TRRparam = [N C0];  
    inital_guess_TRRparam = [bH aH yH bV aV dV EH0 IH0 EV0 IV0]; % variable_TRRparam
    
end

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');

% -----------------------------------------------------------


%% TRR procedure
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- Running TRR Algorithm --- ');
disp(' ');
disp(' ');


% Time of analysis
t0 = 7;                     % initial
t1 = 365;                   % final (starting at january 2016)
dt = 1;                     % time steps (days)

tspan = t0:dt:t1;           % interval of analysis
Ndt   = length(tspan);      % number of time steps


% Output function
    % This function must supply the F(p,xdata) for the ydata comparisson

F = @(variable_TRRparam,day)TRR_FunctionOutput_SEIR_SEI(day,case_name,variable_TRRparam,fixed_TRRparam,tspan,dt,IC); 


% TRR Algorithm

    % Choose bounds
    % Choose time series comparisson. Remember to change on the TRR_FunctionOutput accordingly.
    % result_TRRparam display the new parameters in the order of initial_guess_TRRparam
    
    
    % lsqcurvefit options
    %lb = [1/16.3   1/12   1/8.8   1/11.6     1/10   1/21  0.9*N       0          0        0.99   0  0  ];                % lower bound result of variable parameters
    %up = [1/8      1/3    1/3     1/6.2      1/5    1/11      N       10000      10000    0.999   1  1  ];                % upper bound result of variable parameters
    lb = [0.01  1/15   1/10   0.01     1/12   1/25         0          0           0  0  ];                % lower bound result of variable parameters
    up = [10   1/2    1/2    10      1/3   1/8             N      N      1  1  ];                % upper bound result of variable parameters

    %bH aH yH bV aV dV EH0 IH0 EV0 IV0 (case5)
    
    options = optimoptions('lsqcurvefit','Display','iter','Algorithm','trust-region-reflective','MaxIterations',100000,'MaxFunctionEvaluations',500000,'FunctionTolerance',1e-005,'StepTolerance',1e-005,'FiniteDifferenceType','central');
            % Display: show process/result. 'iter' whill show the steps.
            % Algorithm: choose algorithm. LM is default. LM can't use lb and up
            % MaxIterations: maximum number of iterations. Default is 400
            % MaxFunEvals: maximum number of function evaluations allowed. Default is 100*NumberOfVariables.
    
            
    % Comparing NewCases
    result_TRRparam = lsqcurvefit(F,inital_guess_TRRparam,day,NewCasesData,lb,up,options);
    %disp(['bH = ',num2str(result_TRRparam(1))])
    %disp(['bV = ',num2str(result_TRRparam(4))])
    %disp(['1/aH = ',num2str(1/result_TRRparam(2))])
    %disp(['1/yH = ',num2str(1/result_TRRparam(3))])
    %disp(['1/aV = ',num2str(1/result_TRRparam(5))])
    %disp(['1/dV = ',num2str(1/result_TRRparam(6))])
    %disp(['EH0 = ',num2str(result_TRRparam(7))])
    %disp(['IH0 = ',num2str(result_TRRparam(8))])
    %disp(['EV0 = ',num2str(result_TRRparam(9))])
    %disp(['IV0 = ',num2str(result_TRRparam(10))])
    
toc
% -----------------------------------------------------------


%% Integrate the new initial value problem
% -----------------------------------------------------------
tic
disp(' ');
disp(' --- integration of the nonlinear dynamics with new parameters and IC --- ');
disp(' ');
disp('    ... ');
disp(' ');


% New matrix of parameters
    % New parameter matrix must be the default phys_param except that those corresponding to the variable_param in each 
        % case are replaced by the values of result_param. 
    % Care for the default parameter vector order!
    
% Default TRRparam vector: [N bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 IV0 C0]
% Don't forget to reinstate RH0 between IH0 and SV0 

if strcmp(case_name,'case_config1');
    
    % variable_param = [IV0]
    new_phys_param = [N bH aH yH bV aV dV];  
    new_IC = [SH0 EH0 IH0 RH0 SV0 EV0 result_TRRparam(1) C0];

 elseif strcmp(case_name,'case_config2');

    % variable_param = [EV0 IV0]
    new_phys_param = [N bH aH yH bV aV dV];  
    new_IC = [SH0 EH0 IH0 RH0 SV0 result_TRRparam(1) result_TRRparam(2) C0];

 
elseif strcmp(case_name,'case_config3');

    % variable_TRRparam = [bH bV EV0 IV0]
    new_phys_param = [N result_TRRparam(1) aH yH result_TRRparam(2) aV dV];  
    new_IC = [SH0 EH0 IH0 RH0 SV0 result_TRRparam(3) result_TRRparam(4) C0];


elseif strcmp(case_name,'case_config4');

    % variable_TRRparam = [bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 IV0]
    new_phys_param = [N result_TRRparam(1:6)];
    new_IC = [result_TRRparam(7:9) RH0 result_TRRparam(10:12) C0];

elseif strcmp(case_name,'case_config5');
    
    % variable_TRRparam = [bH aH yH bV aV dV EH0 IH0 EV0 IV0]
    new_phys_param = [N result_TRRparam(1:6)];
    SH0 = N - result_TRRparam(7) - result_TRRparam(8) - 29639;
    SV0 = 1 - result_TRRparam(9) - result_TRRparam(10);
    new_IC = [SH0 result_TRRparam(7:8) RH0 SV0 result_TRRparam(9:10) C0];    
end

% Compartmentalization hypotheses correction
%Xh = N - RH0 - sum(result_TRRparam(7:9))
%Xv = 1 - sum(result_TRRparam(10:12))

%new_IC(1) = result_TRRparam(7) + Xh;
%new_IC(5) = result_TRRparam(10) + Xv;

% Integration with new parameters
        % ODE solver optional parameters
            %opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);
[time y] = ode45(@(t,y)rhs_SEIR_SEI(t,y,new_phys_param),tspan,new_IC);


% Time series 
SH = y(:,1);                % susceptible humans (number of ppl)
EH = y(:,2);                % incubating humans (number of ppl)
IH = y(:,3);                % infectious humans (number of ppl)
RH = y(:,4);                % recovered humans (number of ppl)

SV = y(:,5);                % proportion of susceptible vectors (number of vectors)
EV = y(:,6);                % proportion of incubating vectors (number of vectors)
IV = y(:,7);                % proportion of infectious vectors (number of vectors)

C  = y(:,8);                % cumulative number of infected humans (number of ppl)


% Discrete number of new cases
NewCases = C(1);
for i = 2:52
    NewCases(i) = C( (7*i-7)/dt +1 ) - C( (7*(i-1)-7)/dt +1 );
end

% Sums
SumH = SH(1)+EH(1)+IH(1);

SumV = SV(1)+EV(1)+IV(1);

disp(' '); 
disp([' SumH =  ',num2str(SumH)]);
disp(' '); 
disp(' '); 
disp([' SumV =  ',num2str(SumV)]);
disp(' '); 

% basic reproduction number
R0 = new_phys_param(6)*new_phys_param(2)*new_phys_param(5)/(new_phys_param(7)*new_phys_param(4)*(new_phys_param(6)+new_phys_param(7)))

toc
% -----------------------------------------------------------


%% Plotting
% -----------------------------------------------------------
    % Plots are on a week time scale. Divide 'time' vector by 7. Multiply
        % "day^-1 unit" vectors by 7.
% -----------------------------------------------------------

 
%% Plot number of susceptible humans - SH
% ...........................................................
% 
% gname  = 'SH';                                                          % file name
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
%     set(gca,'position',[0.2 0.2 0.7 0.7]);                              % graph position inside the figure
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
%     gname = [gname, '.eps'];
%     graph_fixPSlinestyle(gname,gname);                                  % fix for grid lines and marker
% 
%     %close;                                                             % close figure
% ...........................................................

%% Plot number of incubating humans - EH
% ...........................................................
% 
% gname  = 'EH';                                                          % file name
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
%     set(gca,'position',[0.2 0.2 0.7 0.7]);                              % graph position inside the figure
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
%     gname = [gname, '.eps'];
%     graph_fixPSlinestyle(gname,gname);                                  % fix for grid lines and marker
%  
%     %close;                                                             % close figure
%...........................................................

%% Plot number of infectious humans - IH
% ...........................................................
% 
% gname  = 'IH';                                                          % file name
% 
% figure('Name',gname,'NumberTitle','off')                                % figure name setting
% 
% fig = plot(time/7,IH,'-b');                                             % plot
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
%                 fs = get(gca,'fontsize');
%                 set(gca,'units','normalized');
%                 xl = xlim;
%                 text(xl(1),yl(2),sprintf('\\times10^{%d}',OM-BD),'fontsize',fs,'VerticalAlignment','bottom');
%     %----------------------------------------------------------------------
%     
%     
%     print(gcf,gname,'-dpdf','-r300','-bestfit');% 
%     %close;                                                             % close figure
%...........................................................

%% Plot number of infectious humans - SumV
% ...........................................................
% 
% gname  = 'SumV';                                                        % file name
% 
% figure('Name',gname,'NumberTitle','off')                                % figure name setting
% 
% fig = plot(time/7,SV+EV+IV,'-b');                                             % plot
%      
%     title('','FontSize',25,'FontName','Helvetica');                     % figure title
%     
%     set(gcf,'color','white');                                           % figure background color
%     
%     xlabel('time (weeks)','FontSize',28,'FontName','Helvetica');        % x axis label and settings
%     ylabel('proportion of vectors','FontSize',28,'FontName','Helvetica');    % y axis label and settings
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
%     set(gca,'FontSize',25,'FontName','Helvetica');                      % tick labels size and font
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
%             BD = 0;                                                     % number of significant figures before the decimal point for the highest tick label
%             OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
%             ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
%             % Formating new tick labels
%                 nyt=cell(size(ryt));
%                 for j=1:length(ryt)
%                     nyt{j}=sprintf('% 1.4f',ryt(j));
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
%     print(gcf,gname,'-dpdf','-r300','-bestfit');% 
%     %close;                                                             % close figure
%...........................................................

%% Plot number of infectious humans - SumH
% ...........................................................
% 
% gname  = 'SumH';                                                        % file name
% 
% figure('Name',gname,'NumberTitle','off')                                % figure name setting
% 
% fig = plot(time/7,SH+EH+IH+RH,'-b');                                    % plot
%      
%     title('','FontSize',25,'FontName','Helvetica');                     % figure title
%     
%     set(gcf,'color','white');                                           % figure background color
%     
%     xlabel('time (weeks)','FontSize',28,'FontName','Helvetica');        % x axis label and settings
%     ylabel('number of people','FontSize',28,'FontName','Helvetica');    % y axis label and settings
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
%     set(gca,'FontSize',25,'FontName','Helvetica');                      % tick labels size and font
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
%             BD = 3;                                                     % number of significant figures before the decimal point for the highest tick label
%             OM = ceil(log10(yl(2)));                                    % ceiling order of magnitude
%             ryt=get(gca,'ytick')/10^(OM-BD);                            % redefining tick labels
%             % Formating new tick labels
%                 nyt=cell(size(ryt));
%                 for j=1:length(ryt)
%                     nyt{j}=sprintf('% 4.1f',ryt(j));
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
%     print(gcf,gname,'-dpdf','-r300','-bestfit');% 
%     %close;                                                             % close figure
%...........................................................

%% Plot cumulative number of infectious humans - C
% ...........................................................

gname  = 'C';                                                           % file name

figure('Name',gname,'NumberTitle','off')                                % figure name setting


fig = plot(time/7,C,'-b');                                              % plot
     
    title('','FontSize',25,'FontName','Helvetica');                     % figure title
    
    set(gcf,'color','white');                                           % figure background color
    
    xlabel('time (weeks)','FontSize',28,'FontName','Helvetica');        % x axis label and settings
    ylabel('cumulative number of people','FontSize',28,'FontName','Helvetica');    % y axis label and settings

    xlim([t0/7 t1/7]);                                                  % x axis limit value
    ylim('auto');                                                       % y axis limit value
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
    
    
    %print(gcf,gname,'-dpdf','-r300','-bestfit');    
    
    %close;                                                             % close figure
% ...........................................................

%% Stem plot new number of cases - NewCases
% ...........................................................

gname  = 'NewCases';                                                    % file name

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
    ylabel('new cases','FontSize',28,'FontName','Helvetica');               % y axis label and settings

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
    
    %print(gcf,gname,'-dpdf','-r300','-bestfit');    
    
    %close 
% ...........................................................     
    
    
