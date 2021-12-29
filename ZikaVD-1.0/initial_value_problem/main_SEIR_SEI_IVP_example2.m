% -----------------------------------------------------------------
% ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
% -----------------------------------------------------------------
% Initial value problem: main_SEIR_SEI_IVP_example2.m
%
% This is the main file for the SEIR-SEI epidemiological model, 
% which a divides the human population in 4 compartments:
%
%   SH = susceptibles humans          (number of individuals)
%   EH = exposed humans               (number of individuals)
%   IH = infectious humans            (number of individuals)
%   RH = recovered humans             (number of individuals)
%
% The vector population is divided in oly 3 compartments:
%
%   SV = proportion of susceptible vectors    (dimensionless)
%   EV = proportion of exposed vectors        (dimensionless)
%   IV = proportion of infectious vectors     (dimensionless)
%
% Infection spreads via cross-transmission between
% a susceptible of one species and a infectious one from the another.
% Delay is modeled as an exposed group: there is an
% latent period until an infected becomes able to 
% transmit (infectious).
% No deaths are considered for humans, all infected become recovered.
% Births and natural deaths are considered for vectors.
%
% The cumulative infectious humans are denoted by:
%
%   C = cumulative number of infectious humans (number of individuals) 
%
% The epidemiological model parameters are:
%
%   N      = human population size (number of individuals)
%
%   betaH  = vector-to-human transmission rate   (days^-1)
%   alphaH = human latent rate                   (days^-1)
%   gamma  = human recovery rate                 (days^-1)
%
%   betaV  = human-to-vector transmission rate   (ind^-1days^-1)
%   alphaV = vector latent rate                  (days^-1)
%   delta  = vector birth/mortality rate         (days^-1)
% 
%
% This code uses rhs_SEIR_SEI.m to define the ODE system
% and outputs the plots and R_nought value. Calculations
% are made on a day time scale.
%
% The conditions simulated by this file was estrated 
% from the Brazil's zika outbreak discrebed in 
% E. Dantas, M. Tosin and A. Cunha Jr,
% Uncertainty quantification in the
% nonlinear dynamics of Zika virus, 
% v3, 2020.
%
% Inputs:
%   param: parameters vector              - double array (7x1)
%   IC: initial conditions vector         - double array (8X1)
%   tspan: time interval                  - double array (?x1)
%   rhs_SEIR_SEI: SEIR-SEI equations file - .m function file
%
% Outputs:
%   R_nought: basic reproduction number   - double
%   figure 1: model state in time         - inplace figure
%   figure 2: number of new cases in time - inplace figure
% -----------------------------------------------------------------
% programmers: Eber Dantas
%              eberdantas@coppe.ufrj.br
%
%              Michel Tosin
%              michel.tosin@uerj.br
%
%              Americo Cunha
%              americo.cunha@uerj.br
%
% number of lines: 114
% last update: Jun 9, 2021
% -----------------------------------------------------------------


clc
clear
close all


% parameters and initial conditions [USER INPUT]
% -----------------------------------------------------------------

% human population size (number of individuals)
N = 206e6; % Brazil population in 2016          

% transmission rates
betaH = 1/8.33; % vector-to-human (days^-1) 
betaV = 1/7.77; % human-to-vector (ind^-1days^-1) 

% latent periods (days)
TalphaH = 12;
TalphaV = 10;
    
% latent rates (days^-1)
alphaH =  1/TalphaH;
alphaV =  1/TalphaV;
        
% human recovery period (days)
Tgamma = 3;

% human recovery rate (days^-1)
gamma  = 1/Tgamma;
    
% Mosquito lifespan (days)
Tdelta = 18;
    
% Mosquito birth/mortality rate (days^-1)    
delta = 1/Tdelta;

% initial conditions
%
% -- Set the initial number of exposed, infectious and recovered humans.
% -- The number of susceptible humans will be the remaining population.
% -- For an invasion scenario, set initial infectious humans to 1.
% -- Set the initial proportion of exposed and infectious mosquitoes.
% -- The proportion of susceptible mosquitoes will be the remaining fraction.

EH0 = 6827;                % initial exposed humans     (number of individuals)
IH0 = 10000;               % initial infectious humans  (number of individuals)
RH0 = 29639;               % initial recovered humans   (number of individuals)
SH0 = N - EH0 - IH0 - RH0; % initial susceptible humans (number of individuals)

EV0 = 4.14e-4;             % initial exposed mosquitoes         (dimensionless)
IV0 = 0;                   % initial infectious mosquitoes      (dimensionless)
SV0 = 1 - EV0 - IV0;       % initial susceptible mosquitoes     (dimensionless)

C0 = 8201; % initial cumulative infectious humans (number of individuals)



% display program header on screen
% -----------------------------------------------------------------

% computing the basic reproduction number R_nought
R_nought = (betaH*alphaV*betaV)/(gamma*(alphaV+delta)*delta);

disp('                                                            ')
disp('============================================================')
disp('             ARBO - Arbovirus Modeling and                  ')
disp('          Uncertainty Quantification Toolbox                ')
disp('                                                            ')
disp('            by Eber Dantas, Michel Tosin,                   ')
disp('       Americo Cunha Jr and Rebecca E. Morrison             ')
disp('                                                            ')
disp('     This is a package for simulation and analysis          ')
disp('            of arbovirus nonlinear dynamics.                ')
disp('============================================================')
disp('                                                            ')
disp(' -----------------------------------------------------------')
disp(' +++++++++++++++++++++ SEIR-SEI model +++++++++++++++++++++ ')
disp(' -----------------------------------------------------------')
disp(['  * human total population =            ',num2str(N)]       )
disp( '    (individuals)       '                                   )
disp(['  * vector-to-human transmission rate = ',num2str(betaH)]   )
disp( '    (days^-1)           '                                   )
disp(['  * human latent rate                 = ',num2str(alphaH)]  )
disp( '    (days^-1)           '                                   )
disp(['  * human recovery rate               = ',num2str(gamma)]   )
disp( '    (days^-1)           '                                   )
disp(['  * human-to-vector transmission rate = ',num2str(betaV)]   )
disp( '    (ind^-1days^-1)           '                             )
disp(['  * vector latent rate                = ',num2str(alphaH)]  )
disp( '    (days^-1)           '                                   )
disp(['  * vector birth/mortality rate       = ',num2str(delta)]   )
disp( '    (days^-1)           '                                   )
disp(['  * R_nought                          = ',num2str(R_nought)])
disp( '    (dimensionless)     '                                   )
disp(' -----------------------------------------------------------')
% -----------------------------------------------------------------


% integration of the initial value problem
% -----------------------------------------------------------------

% parameters vector
param = [N betaH alphaH gamma betaV alphaV delta];

% initial conditions vector
IC = [SH0 EH0 IH0 RH0 SV0 EV0 IV0 C0];

% time interval of analysis
   t0 = 7;                  % initial time (days)
   t1 = 364;                % final time   (days)
   dt = 1  ;                % time steps   (days)
tspan = t0:dt:t1;           % interval of analysis


% ODE solver Runge-Kutta45
[time, y] = ode45(@(t,y)rhs_SEIR_SEI(t,y,param),tspan,IC);

% Time series 
SH = y(:,1); % susceptible humans                     (number of ind)
EH = y(:,2); % exposed humans                         (number of ind)
IH = y(:,3); % infectious humans                      (number of ind)
RH = y(:,4); % recovered humans                       (number of ind)

SV = y(:,5); % proportion of susceptible vectors      (dimensionless)
EV = y(:,6); % proportion of exposed vectors          (dimensionless)
IV = y(:,7); % proportion of infectious vectors       (dimensionless)

C  = y(:,8); % cumulative number of infectious humans (number of ind)
% -----------------------------------------------------------------

% post-processing
% -----------------------------------------------------------------

% NewCases (per week) computation
tu = 7;                     % time unit (days)
%
% -- tu/dt must be an integer
% -- tu = 1 defines a 'per day' computation
% -- For a 'per week' computation, use tu = 7 and t0 = 7

h = round(tu/dt); % step of convertion
Ctu = C(1:h:end); % extract the C in the tu instants
NC = [Ctu(1); (Ctu(2:end)-Ctu(1:end-1) )]; % New number of infectious humans

% custom color
yellow = [255 204  0]/256;

% plot all humans compartments
figure(1)
hold on
fig1(1) = plot(time,SH);
fig1(2) = plot(time,EH);
fig1(3) = plot(time,IH);
fig1(4) = plot(time,RH);
fig1(5) = plot(time,C);
hold off

    % plot labels
     title('SEIR-SEI human compartments');
    xlabel('time (days)'                );
    ylabel('number of individuals'      );

    % set plot settings
    set(gca,'FontSize',18);
    set(fig1,{'Color'},{'b';yellow;'r';'g';'m'});
    set(fig1,{'LineWidth'},{3;3;3;3;3});

    % legend
    leg = {'Suceptibles'; 'Exposed'; 'Infectious'; 'Recovered'; 'Cum. Infectious'};
    legend(fig1,leg,'FontSize',15,'location','northeast');
    
    % axis limits
    xlim([t0 t1]);
    ylim([0 10000]);

    saveas(figure(1),'fig_IVP_example2_human_compart.png')

% plot all humans compartments
figure(2)
hold on
fig2(1) = plot(time,SV);
fig2(2) = plot(time,EV);
fig2(3) = plot(time,IV);
hold off

    % plot labels
     title('SEIR-SEI vectors compartments');
    xlabel('time (days)'                  );
    ylabel('proportion of vectors'        );

    % set plot settings
    set(gca,'FontSize',18);
    set(fig2,{'Color'},{'b';yellow;'r'});
    set(fig2,{'LineWidth'},{3;3;3});

    % legend
    leg = {'Suceptibles'; 'Exposed'; 'Infectious'};
    legend(fig2,leg,'FontSize',15,'location','northeast');
    
    % axis limits
    xlim([t0 t1]);
    ylim([0 1e-3]);

    saveas(figure(2),'fig_IVP_example2_vector_compart.png')

% plot NewCases (per week) of SEIR-SEI model
figure(3)
fig3 = scatter([1:length(NC)],NC);

    % plot labels
     title('New Cases per week (SEIR-SEI)');
    xlabel('time (weeks)'                 );
    ylabel('number of individuals'        );

    % set plot settings
    set(gca,'FontSize',18);
    set(fig3,{'MarkerFaceColor','MarkerEdgeColor'},{'y','r'});

    % axis limits
    xlim([1 length(NC)]);
    
    saveas(figure(3),'fig_IVP_example2_NC.png')
% -----------------------------------------------------------------
