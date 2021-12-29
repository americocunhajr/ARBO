% -----------------------------------------------------------------
% ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
% -----------------------------------------------------------------
% Uncertainty quantification: main_SEIR_SEI_MC_example1.m
%
% This is the main file for a uncertainty propagation simulation 
% for the SEIR-SEI epidemiological model. Details about the model 
% can be found in the Initial value problem module.
%
% The random input is the parameter
%
%   delta  = vector birth/mortality rate         (days^-1)
%
% considering the input information
%
%   delta_min = minimum value for the parameter
%   delta_max = maximum value for the parameter
%
% This code uses rhs_SEIR_SEI.m to define the ODE system
% and outputs the plots of confidente bands for several model outputs. 
% Calculations are made on a day time scale.
%
% The nominal conditions simulated by this file was estrated 
% from the Brazil's zika outbreak discrebed in 
% E. Dantas, M. Tosin and A. Cunha Jr,
% Calibration of a SEIR–SEI epidemic model 
% to describe the Zika virus outbreak in Brazil, 
% Applied Mathematics and Computation, 338, p.249-259, 
% 2020. DOI: https://doi.org/10.1016/j.amc.2018.06.024
%
% Inputs:
%   Ns: Number of samples                 - double
%   param: parameters vector              - double array (7x1)
%   IC: initial conditions vector         - double array (8X1)
%   tspan: time interval                  - double array (?x1)
%   rhs_SEIR_SEI: SEIR-SEI equations file - .m function file
%
% Outputs:
%   figure 1: MC convergence metric       - inplace figure
%   figure 2: SH confidence band          - inplace figure
%   figure 3: EH confidence band          - inplace figure
%   figure 4: IH confidence band          - inplace figure
%   figure 5: RH confidence band          - inplace figure
%   figure 6: C confidence band           - inplace figure
%   figure 7: NC confidence band          - inplace figure
% -----------------------------------------------------------------
% programmers: Michel Tosin
%              michel.tosin@uerj.br
%
%              Eber Dantas
%              eberdantas@coppe.ufrj.br
%
%              Americo Cunha
%              americo.cunha@uerj.br
%
% number of lines: 228
% last update: Jun 23, 2021
% -----------------------------------------------------------------


clc
clear
close all


% display program header on screen
% -----------------------------------------------------------------
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
% -----------------------------------------------------------------


%% ----------------------------------------------------------------
% Pre-processing
%% ----------------------------------------------------------------


% Data - Brazil
% -----------------------------------------------------------------
% Load 2016 data
ProbData = load('2016zika_Prob250117.dat');

% Data time vector
day = ProbData(1,:)*7;

% New cases data
NCData = ProbData(2,:);

% Cumulative cases data
CData = cumsum(NCData);

% Load 2015 data - Cumulative number of cases
C2015 = 29639;      % until EW 49
% -----------------------------------------------------------------


% parameters and initial conditions [USER INPUT]
% -----------------------------------------------------------------
% human population size (number of individuals)
N = 206e6; % Brazil population in 2016          

% transmission rates
betaH = 1/8.33; % vector-to-human (days^-1) 
betaV = 1/7.77;  % human-to-vector (ind^-1days^-1) 

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
C0 = CData(1); % initial cumulative infectious humans (number of individuals)

EH0 = 6827;                % initial exposed humans     (number of individuals)
IH0 = 10000;               % initial infectious humans  (number of individuals)
RH0 = C2015;               % initial recovered humans   (number of individuals)
SH0 = N - EH0 - IH0 - RH0; % initial susceptible humans (number of individuals)

EV0 = 4.14e-4;             % initial exposed mosquitoes         (dimensionless)
IV0 = 0;                   % initial infectious mosquitoes      (dimensionless)
SV0 = 1 - EV0 - IV0;       % initial susceptible mosquitoes     (dimensionless)
% -----------------------------------------------------------------


% nominal model
% -----------------------------------------------------------------
% time interval of analysis
   t0 = 7;                  % initial time (days)
   t1 = 364;                % final time   (days)
   dt = 1;                  % time steps   (days)
tspan = t0:dt:t1;           % interval of analysis
  Ndt = length(tspan);      % number of time instants
   tu = 7;                  % unit for weeks
    h = round(tu/dt);       % step for unit convertion  
   
    % parameters vector
    param = [N betaH alphaH gamma betaV alphaV delta];

    % initial conditions
    IC = [SH0 EH0 IH0 RH0 SV0 EV0 IV0 C0];

    % ODE solver Runge-Kutta45
    [time,y] = ode45(@(t,y)rhs_SEIR_SEI(t,y,param),tspan,IC);

    % time series of susceptible humans (number of ppl)
    SH_Nominal = y(:,1)';

    % time series of incubating humans (number of ppl)
    EH_Nominal = y(:,2)';

    % time series of infectious humans (number of ppl)
    IH_Nominal = y(:,3)';

    % time series of recovered humans (number of ppl)
    RH_Nominal = y(:,4)';

    % time series of proportion of susceptible vectors (number of vectors)
    SV_Nominal = y(:,5)';

    % time series of proportion of incubating vectors (number of vectors)
    EV_Nominal = y(:,6)';

    % time series of proportion of infectious vectors (number of vectors)
    IV_Nominal = y(:,7)';

    % time series of cumulative number of infected humans (number of ppl)
     C_Nominal = y(:,8)';
    
    % discrete number of new cases
    Ctu = C_Nominal(1:h:end);       % extract the C in the tu instants
    NC_Nominal = [Ctu(1) ( Ctu(2:end) - Ctu(1:end-1) )]; 
    
% -----------------------------------------------------------------


% define stochastic parameters
% -----------------------------------------------------------------
% RNG seed
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream); % Matlab 2013

% number of random samples
%Ns = 2;
%Ns = 4;
%Ns = 16;
%Ns = 64;
%Ns = 256;
Ns = 1024;
%Ns = 4096;
%Ns = 16384;
%Ns = 65536;

% known information for param
delta_min = 1/21;   % min value
delta_max = 1/11;   % max value

% generate samples for betaH and delta
delta_samp = delta_min + (delta_max-delta_min).*rand(Ns,1);
% -----------------------------------------------------------------


%% ----------------------------------------------------------------
% Processing
%% ----------------------------------------------------------------


% Monte Carlo simulation
% -----------------------------------------------------------------
% preallocate memory for system QoIs
SH_MC = zeros(Ns,Ndt);
EH_MC = zeros(Ns,Ndt);
IH_MC = zeros(Ns,Ndt);
RH_MC = zeros(Ns,Ndt);
SV_MC = zeros(Ns,Ndt);
EV_MC = zeros(Ns,Ndt);
IV_MC = zeros(Ns,Ndt);
 C_MC = zeros(Ns,Ndt);

NC_MC = zeros(Ns,length(NC_Nominal));

% preallocate memory for MC conv metric
MC_conv = zeros(1,Ns);

for imc=1:Ns
    
    if mod(imc,sqrt(Ns)) == 0
        disp('')
        disp(imc)
    end

    % imc-th param vector
    param(7) = delta_samp(imc);

    % ODE solver Runge-Kutta45
    [time,y] = ode45(@(t,y)rhs_SEIR_SEI(t,y,param),tspan,IC);

    % time series of susceptible humans (number of ppl)
    SH_MC(imc,:) = y(:,1);

    % time series of incubating humans (number of ppl)
    EH_MC(imc,:) = y(:,2);

    % time series of infectious humans (number of ppl)
    IH_MC(imc,:) = y(:,3);

    % time series of recovered humans (number of ppl)
    RH_MC(imc,:) = y(:,4);

    % time series of proportion of susceptible vectors (number of vectors)
    SV_MC(imc,:) = y(:,5);

    % time series of proportion of incubating vectors (number of vectors)
    EV_MC(imc,:) = y(:,6);

    % time series of proportion of infectious vectors (number of vectors)
    IV_MC(imc,:) = y(:,7);

    % time series of cumulative number of infected humans (number of ppl)
     C_MC(imc,:) = y(:,8);
    
    % compute MC convergence metric
    MC_conv(1,imc) = trapz(time,SH_MC(imc,:).^2 + EH_MC(imc,:).^2 + ...
                                IH_MC(imc,:).^2 + RH_MC(imc,:).^2 + ...
                                SV_MC(imc,:).^2 + EV_MC(imc,:).^2 + ...
                                IV_MC(imc,:).^2 +  C_MC(imc,:).^2);

end

% compute MC convergence metric
MC_conv = randvar_mc_conv(MC_conv,Ns);

% discrete number of new cases
Ctu = C_MC(:,1:h:end);       % extract the C in the tu instants
NC_MC = [Ctu(:,1) ( Ctu(:,2:end) - Ctu(:,1:end-1) )];
% -----------------------------------------------------------------


%% ----------------------------------------------------------------
% Post-processing
%% ----------------------------------------------------------------


% compute the statistics
% -----------------------------------------------------------------
% sample average of random processes
SH_smp_avg = mean(SH_MC);
EH_smp_avg = mean(EH_MC);
IH_smp_avg = mean(IH_MC);
RH_smp_avg = mean(RH_MC);
SV_smp_avg = mean(SV_MC);
EV_smp_avg = mean(EV_MC);
IV_smp_avg = mean(IV_MC);
 C_smp_avg = mean(C_MC);
NC_smp_avg = mean(NC_MC);

% confidence probability (percentual)
Pc = 95;

% upper percentil
r_plus = 0.5*(100 + Pc);

% lower percentil
r_minus = 0.5*(100 - Pc);

% confidence bands upper bounds
SH_upp = prctile(SH_MC,r_plus);
EH_upp = prctile(EH_MC,r_plus);
IH_upp = prctile(IH_MC,r_plus);
RH_upp = prctile(RH_MC,r_plus);
SV_upp = prctile(SV_MC,r_plus);
EV_upp = prctile(EV_MC,r_plus);
IV_upp = prctile(IV_MC,r_plus);
 C_upp = prctile( C_MC,r_plus);
NC_upp = prctile(NC_MC,r_plus);

% confidence bands lower bounds
SH_low = prctile(SH_MC,r_minus);
EH_low = prctile(EH_MC,r_minus);
IH_low = prctile(IH_MC,r_minus);
RH_low = prctile(RH_MC,r_minus);
SV_low = prctile(SV_MC,r_minus);
EV_low = prctile(EV_MC,r_minus);
IV_low = prctile(IV_MC,r_minus);
 C_low = prctile( C_MC,r_minus);
NC_low = prctile(NC_MC,r_minus);   
% -----------------------------------------------------------------

% Plots
% -----------------------------------------------------------------
% plot MC convergence
figure(1)
hold on
fig1 = plot(1:Ns,MC_conv,'xb');

    % plot labels
     title('Monte Carlo convergence'  );
    xlabel('number of MC realizations');
    ylabel('convergence metric'       );

    % set plot settings
    set(gca,'FontSize',18);
    set(fig1,'LineWidth',3);

    % axis limits
    xlim([1 Ns]);
    ylim('auto');
    
    saveas(figure(1),'fig_UQ_example1_MCconv.png')

% plot SH confidence band
figure(2)
hold on
fig2(1) = patch([time' fliplr(time')],[SH_upp fliplr(SH_low)],[0.85 0.85 0.85],'EdgeColor','none');
fig2(2) = plot(time,SH_smp_avg,'b'); 
fig2(3) = plot(time,SH_Nominal,'-g'); 
hold off

    % plot labels
     title('SH confidence band'   );
    xlabel('time (days)'          );
    ylabel('number of individuals');

    % set plot settings
    set(gca,'FontSize',18);
    set(fig2,{'LineWidth'},{3;3;3});

    % legend
    leg = {[num2str(Pc),'% prob.']; 'mean'; 'nominal'};
    legend(fig2,leg,'FontSize',15,'location','northeast');
    
    % axis limits
    xlim([day(1) day(end)]);
    ylim('auto');
    
    saveas(figure(2),'fig_UQ_example1_SHband.png')

% plot EH confidence band
figure(3)
hold on
fig3(1) = patch([time' fliplr(time')],[EH_upp fliplr(EH_low)],[0.85 0.85 0.85],'EdgeColor','none');
fig3(2) = plot(time,EH_smp_avg,'b'); 
fig3(3) = plot(time,EH_Nominal,'-g'); 
hold off

    % plot labels
     title('EH confidence band'   );
    xlabel('time (days)'          );
    ylabel('number of individuals');

    % set plot settings
    set(gca,'FontSize',18);
    set(fig3,{'LineWidth'},{3;3;3});

    % legend
    leg = {[num2str(Pc),'% prob.']; 'mean'; 'nominal'};
    legend(fig3,leg,'FontSize',15,'location','northeast');
    
    % axis limits
    xlim([day(1) day(end)]);
    ylim('auto');
    
    saveas(figure(3),'fig_UQ_example1_EHband.png')

% plot IH confidence band
figure(4)
hold on
fig4(1) = patch([time' fliplr(time')],[IH_upp fliplr(IH_low)],[0.85 0.85 0.85],'EdgeColor','none');
fig4(2) = plot(time,IH_smp_avg,'b'); 
fig4(3) = plot(time,IH_Nominal,'-g'); 
hold off

    % plot labels
     title('IH confidence band'   );
    xlabel('time (days)'          );
    ylabel('number of individuals');

    % set plot settings
    set(gca,'FontSize',18);
    set(fig4,{'LineWidth'},{3;3;3});

    % legend
    leg = {[num2str(Pc),'% prob.']; 'mean'; 'nominal'};
    legend(fig4,leg,'FontSize',15,'location','northeast');
    
    % axis limits
    xlim([day(1) day(end)]);
    ylim('auto');
    
    saveas(figure(4),'fig_UQ_example1_IHband.png')

% plot RH confidence band
figure(5)
hold on
fig5(1) = patch([time' fliplr(time')],[RH_upp fliplr(RH_low)],[0.85 0.85 0.85],'EdgeColor','none');
fig5(2) = plot(time,RH_smp_avg,'b'); 
fig5(3) = plot(time,RH_Nominal,'-g'); 
hold off

    % plot labels
     title('RH confidence band'   );
    xlabel('time (days)'          );
    ylabel('number of individuals');

    % set plot settings
    set(gca,'FontSize',18);
    set(fig5,{'LineWidth'},{3;3;3});

    % legend
    leg = {[num2str(Pc),'% prob.']; 'mean'; 'nominal'};
    legend(fig5,leg,'FontSize',15,'location','northeast');
    
    % axis limits
    xlim([day(1) day(end)]);
    ylim('auto');
    
    saveas(figure(5),'fig_UQ_example1_RHband.png')

% plot C confidence band
figure(6)
hold on
fig6(1) = patch([time' fliplr(time')],[C_upp fliplr(C_low)],[0.85 0.85 0.85],'EdgeColor','none');
fig6(2) = plot(time,C_smp_avg,'b'); 
fig6(3) = plot(time,C_Nominal,'-g'); 
fig6(4) = scatter(day,CData,'m');
hold off

    % plot labels
     title('C confidence band'    );
    xlabel('time (days)'          );
    ylabel('number of individuals');

    % set plot settings
    set(gca,'FontSize',18);
    set(fig6,{'LineWidth'},{3;3;3;3});

    % legend
    leg = {[num2str(Pc),'% prob.']; 'mean'; 'nominal'; 'data'};
    legend(fig6,leg,'FontSize',15,'location','northeast');
    
    % axis limits
    xlim([day(1) day(end)]);
    ylim('auto');
    
    saveas(figure(6),'fig_UQ_example1_Cband.png')
    
% plot NC confidence band
figure(7)
hold on
fig7(1) = patch([day/7 fliplr(day/7)],[NC_upp fliplr(NC_low)],[0.85 0.85 0.85],'EdgeColor','none');
fig7(2) = plot(day/7,NC_smp_avg,'b'); 
fig7(3) = plot(day/7,NC_Nominal,'-g'); 
fig7(4) = scatter(day/7,NCData,'m');
hold off

    % plot labels
     title('NC confidence band'   );
    xlabel('time (weeks)'         );
    ylabel('number of individuals');

    % set plot settings
    set(gca,'FontSize',18);
    set(fig7,{'LineWidth'},{3;3;3;3});

    % legend
    leg = {[num2str(Pc),'% prob.']; 'mean'; 'nominal'; 'data'};
    legend(fig7,leg,'FontSize',15,'location','northeast');
    
    % axis limits
    xlim([day(1)/7 day(end)/7]);
    ylim('auto');  

    saveas(figure(7),'fig_UQ_example1_NCband.png')
% -----------------------------------------------------------------
