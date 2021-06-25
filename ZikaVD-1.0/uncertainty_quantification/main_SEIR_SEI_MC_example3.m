% -----------------------------------------------------------------
% ZikaVD - Zika Virus Dynamics
% -----------------------------------------------------------------
% Uncertainty quantification: main_SEIR_SEI_MC_example2.m
%
% This is the main file for a uncertainty propagation simulation 
% for the SEIR-SEI epidemiological model. Details about the model 
% can be found in the Initial value problem module.
%
% The random input is the Initial condition
%
%   IH0  = initial number of infectious humans  (number of indiv)
%
% considering the input information
%
%   IH0_min = minimum value for the IC
%   IH0_max = maximum value for the IC
%   IH0_mu = mean value for the IC
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
%   figure 2: C confidence band           - inplace figure
%   figure 3: NC confidence band          - inplace figure
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
% number of lines: 172
% last update: Jun 23, 2021
% -----------------------------------------------------------------


clc
clear
close all


% display program header on screen
% -----------------------------------------------------------------
disp('                                                            ')
disp('============================================================')
disp('   ZikaVD - Zika Virus Dynamics                             ')
disp('   by Michel Tosin, Eber Dantas and Americo Cunha Jr        ')
disp('                                                            ')
disp('   This is an easy-to-run Matlab code to simulate           ')
disp('   the nonlinear dynamics of the Zika virus.                ')
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

% known information for the param
%            betaH  delta
param_min = [1/16.3 1/21];   % min value
param_max = [ 1/8   1/11];   % max value
param_mu  = [betaH delta];   % mean value (1st moment)

% uniform samples in (0,1)
U = rand(Ns,2);

% number of mesh points for random variable support
Nx = 1000;

% Lagrange multipliers for MaxEnt distribution
for i = 1:2
    [lambda(:,i),pdf_param(:,i),supp_param(:,i)] = ...
      maxent_lagrange_mc(param_min(i),param_max(i),Nx,[1; param_mu(i)]);
    l0 = lambda(1,i);
    l1 = lambda(2,i);

    % CDF of the i-th param
    cdf = @(x) exp(-l0)*(exp(-l1*param_min(i))-exp(-l1*x))/l1;

    % quantile function of Alpha
    cdf_inv = @(x) -log(exp(-l1*param_min(i))-l1*x*exp(l0))/l1;

    % MaxEnt random samples for the i-th param
    param_samp(:,i) = cdf_inv(U(:,i));
    param_samp(:,i) = real(param_samp(:,i));
end
% -----------------------------------------------------------------


%% ----------------------------------------------------------------
% Processing
%% ----------------------------------------------------------------


% Monte Carlo simulation
% -----------------------------------------------------------------
% preallocate memory for system QoIs
 C_MC = zeros(Ns,Ndt);
NC_MC = zeros(Ns,length(NC_Nominal));

% preallocate memory for MC conv metric
MC_conv = zeros(1,Ns);

for imc=1:Ns
    
    if mod(imc,sqrt(Ns)) == 0
        disp('')
        disp(imc)
    end

    % imc-th IC vector
    param(2) = param_samp(imc,1);
    param(7) = param_samp(imc,2);

    % ODE solver Runge-Kutta45
    [time,y] = ode45(@(t,y)rhs_SEIR_SEI(t,y,param),tspan,IC);

    % time series of cumulative number of infected humans (number of ppl)
     C_MC(imc,:) = y(:,8);
    
    % compute MC convergence metric
    MC_conv(1,imc) = trapz(time,C_MC(imc,:).^2);

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
 C_smp_avg = mean(C_MC);
NC_smp_avg = mean(NC_MC);

% time average of random processes
 C_time_avg = mean(C_MC,2);
NC_time_avg = mean(NC_MC,2);

% standard deviation of random processes
 C_std = std(C_MC);
NC_std = std(NC_MC);

% skewness of random processes
 C_skew = skewness(C_MC);
NC_skew = skewness(NC_MC);

% kurtosis of the random processes
 C_kurt = kurtosis(C_MC)-3;
NC_kurt = kurtosis(NC_MC)-3;

% confidence probability (percentual)
Pc = 95;

% upper percentil
r_plus = 0.5*(100 + Pc);

% lower percentil
r_minus = 0.5*(100 - Pc);

% confidence bands upper bounds
 C_upp = prctile( C_MC,r_plus);
NC_upp = prctile(NC_MC,r_plus);

% confidence bands lower bounds
 C_low = prctile( C_MC,r_minus);
NC_low = prctile(NC_MC,r_minus);

% number of bins
Nbins = round(sqrt(Ns));

% number of KSD points
Nksd = 100;

% time-average histograms estimation
[C_bins,  C_freq] = randvar_pdf(randvar_normalize(C_time_avg),Nbins);
[NC_bins,NC_freq] = randvar_pdf(randvar_normalize(NC_time_avg),Nbins);

% time-average PDFs curves estimation
[C_ksd,  C_supp] = randvar_ksd(randvar_normalize(C_time_avg),Nksd);
[NC_ksd,NC_supp] = randvar_ksd(randvar_normalize(NC_time_avg),Nksd);

% random variables statistics
C_time_avg_mean =     mean(C_time_avg)
C_time_avg_std  =      std(C_time_avg)
C_time_avg_skew = skewness(C_time_avg)
C_time_avg_kurt = kurtosis(C_time_avg)-3

C_time_avg_low = prctile(C_time_avg,r_minus)
C_time_avg_upp = prctile(C_time_avg,r_plus)

NC_time_avg_mean =     mean(NC_time_avg)
NC_time_avg_std  =      std(NC_time_avg)
NC_time_avg_skew = skewness(NC_time_avg)
NC_time_avg_kurt = kurtosis(NC_time_avg)-3

NC_time_avg_low = prctile(NC_time_avg,r_minus)
NC_time_avg_upp = prctile(NC_time_avg,r_plus)

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

    saveas(figure(1),'fig_UQ_example3_MCconv.png')

% plot C confidence band
figure(2)
hold on
fig2(1) = patch([time' fliplr(time')],[C_upp fliplr(C_low)],[0.85 0.85 0.85],'EdgeColor','none');
fig2(2) = plot(time,C_smp_avg,'b'); 
fig2(3) = plot(time,C_Nominal,'-g'); 
fig2(4) = scatter(day,CData,'m');
hold off

    % plot labels
     title('C confidence band'    );
    xlabel('time (days)'          );
    ylabel('number of individuals');

    % set plot settings
    set(gca,'FontSize',18);
    set(fig2,{'LineWidth'},{3;3;3;3});

    % legend
    leg = {[num2str(Pc),'% prob.']; 'mean'; 'nominal'; 'data'};
    legend(fig2,leg,'FontSize',15,'location','northeast');
    
    % axis limits
    xlim([day(1) day(end)]);
    ylim('auto');
    
    saveas(figure(2),'fig_UQ_example3_Cband.png')
    
% plot NC confidence band
figure(3)
hold on
fig3(1) = patch([day/7 fliplr(day/7)],[NC_upp fliplr(NC_low)],[0.85 0.85 0.85],'EdgeColor','none');
fig3(2) = plot(day/7,NC_smp_avg,'b'); 
fig3(3) = plot(day/7,NC_Nominal,'-g'); 
fig3(4) = scatter(day/7,NCData,'m');
hold off

    % plot labels
     title('NC confidence band'   );
    xlabel('time (weeks)'         );
    ylabel('number of individuals');

    % set plot settings
    set(gca,'FontSize',18);
    set(fig3,{'LineWidth'},{3;3;3;3});

    % legend
    leg = {[num2str(Pc),'% prob.']; 'mean'; 'nominal'; 'data'};
    legend(fig3,leg,'FontSize',15,'location','northeast');
    
    % axis limits
    xlim([day(1)/7 day(end)/7]);
    ylim('auto');      

    saveas(figure(3),'fig_UQ_example3_NCband.png')
% -----------------------------------------------------------------
