
% -----------------------------------------------------------
%  main_zikav_SEIR_SEI_mc.m
%
%  This script is the main file for a program that simulates
%  the stochastic nonlinear dynamics of Zika virus using a 
%  SEIR-SEI epidemic model with 7 equations and 8 parameters.
%
% Calculations are made on a day time scale and plotting 
% on a week time scale
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
%               americo.cunhajr@gmail.com
%               
%  last update: March 1, 2019
% -----------------------------------------------------------------


clc
clear all
close all


% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' ---------------------------------------------------')
disp(' Zika Virus Nonlinear Dynamics (SEIR-SEI model)     ')
disp(' (uncertanty propagation calculation)               ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Barbosa da Cunha Junior                    ')
disp(' americo.cunhajr@gmail.com                          ')
disp(' ---------------------------------------------------')
disp('                                                    ')
% -----------------------------------------------------------


% simulation information
% -----------------------------------------------------------
case_name = 'zikav_SEIR_SEI_mc';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------


% Data - Brazil
% -----------------------------------------------------------------
disp(' '); 
disp(' --- load data --- ');
disp(' ');
disp('    ... ');
disp(' ');

% Load 2016 data
ProbData = load('2016zika_Prob250117.dat');
ConfData = load('2016zika_Conf250117.dat');

% Data time vector
day = ProbData(1,:)*7;
%day = ConfData(1,:)*7;

% New cases data
NewCasesData = ProbData(2,:);
%NewCasesData = ConfData(2,:);

% Cumulative cases data
CData = cumsum(NewCasesData);

% Load 2015 data - Cumulative number of cases
%C2015 = 27127;     % until EW 45
C2015 = 29639;      % until EW 49

%C2015 = 443502;     % inferior PAHO estimation     
%C2015 = 1301140;    % superior PAHO estimation
% -----------------------------------------------------------


% define deterministic physical parameters
% -----------------------------------------------------------
disp(' '); 
disp(' --- define model parameters --- ');
disp(' ');
disp('    ... ');
disp(' ');


% Human population size (number of ppl)
N = 206.0e6;

% Intrinsic incubation period (days)
TaH = 12.0;

% Extrinsic incubation period (days)
TaV = 10.0;
        
% Human infectious period (days)
TyH = 3.0;
        
% Mosquito lifespan (days)
TdV = 21.0;

% Rates and inverses
aH  = 1./TaH;  % intrinsic incubation ratio (days^-1)
aV  = 1./TaV;  % extrinsic incubation ratio (days^-1)
yH  = 1./TyH;  % human infectious ratio (days^-1)
dV  = 1./TdV;  % mosquito death rate (lifespan inverse) (days^-1)
    
% Transmission rates
bH  = 1/10.40;   % vector-to-human transmission rate (days^-1)
bV  = 1/7.77 ;   % human-to-vector transmission rate (days^-1)

% Initial conditions
IH0 = 10000;                  % initial number infectious humans
EH0 = 6827;                   % initial number incubating humans   
SH0 = N - IH0 - EH0 - C2015;  % initial number susceptible humans
RH0 = C2015;                  % initial number recovered humans        

IV0 = 0;                      % initial proportion infectious vectors
EV0 = 0.000414;               % initial proportion incubating vectors
SV0 = 1.0 - IV0 - EV0;        % initial proportion susceptible vectors

C0 = NewCasesData(1);         % initial cumulative infectious humans
% -----------------------------------------------------------


% nominal model
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- integrate nominal model --- ');
disp(' ');
disp('    ... ');
disp(' ');

    % initial time of analysis
    t0 = 7;

    % final time of analysis
    % (analysis starting at january 2016)
    t1 = 365;
    
    % value of time steps
    dt = 1;

    % ODE solver optional parameters
    %opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);

    % interval of analysis
    tspan = t0:dt:t1;

    % number of time steps
    Ndt = length(tspan);

    % physical parameters vector
    phys_param = [N bH aH yH bV aV dV];

    % initial conditions
    IC = [SH0 EH0 IH0 RH0 SV0 EV0 IV0 C0];

    % ODE solver Runge-Kutta45
    [time,y] = ode45(@(t,y)rhs_SEIR_SEI(t,y,phys_param),tspan,IC);

    % time series of susceptible humans (number of ppl)
    Nominal_SH = y(:,1);

    % time series of incubating humans (number of ppl)
    Nominal_EH = y(:,2);

    % time series of infectious humans (number of ppl)
    Nominal_IH = y(:,3);

    % time series of recovered humans (number of ppl)
    Nominal_RH = y(:,4);

    % time series of proportion of susceptible vectors (number of vectors)
    Nominal_SV = y(:,5);

    % time series of proportion of incubating vectors (number of vectors)
    Nominal_EV = y(:,6);

    % time series of proportion of infectious vectors (number of vectors)
    Nominal_IV = y(:,7);

    % time series of cumulative number of infected humans (number of ppl)
    Nominal_C = y(:,8);
    
    % discrete number of new cases
    Nweeks = 52;

    Nominal_NewCases = zeros(1,Nweeks);

    Nominal_NewCases(1,1) = C0;
    
    for inw = 2:Nweeks
        Nominal_NewCases(1,inw) = ...
            Nominal_C((7*inw-7)/dt+1) - Nominal_C((7*(inw-1)-7)/dt+1);
    end
    
toc
% -----------------------------------------------------------


% define stochastic parameters
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- generate random samplss --- ');
disp(' ');
disp('    ... ');
disp(' ');

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

% known information for IH0
IH0_min   = 10.0e3;
IH0_max   = 25.0e3;
IH0_mu    = IH0;
IH0_delta = 0.25;
IH0_sig   = IH0_mu*IH0_delta;
IH0_mu2   = IH0_mu^2 + IH0_sig^2;

% number of mesh points for random variable support
Nx = 1000;

% Lagrange multipliers for MaxEnt distribution
[lambda,pdf_IH0,supp_IH0] = ...
    maxent_lagrange_mc(IH0_min,IH0_max,Nx,[1; IH0_mu; IH0_mu2]);
l0 = lambda(1);
l1 = lambda(2);
l3 = lambda(3);

% CDF of IH0
cdf_IH0 = @(x) exp(-l0)*(exp(-l1*IH0_min)-exp(-l1*x))/l1;

% quantile function of Alpha
cdf_inv_IH0 = @(x) -log(exp(-l1*IH0_min)-l1*x*exp(l0))/l1;

% uniform samples in (s0,1)
U = rand(Ns,1);

% MaxEnt random samples for IH0 (inverse transform method)
IH0_samples = cdf_inv_IH0(U);
IH0_samples = real(IH0_samples);

% generate samples for human initial conditions
%IH0 = IH0min + (IH0max-IH0min)*rand(Ns,1);
%EH0 = rand(Ns,1).*(N-IH0);
%RH0 = rand(Ns,1).*(N-IH0-EH0);

toc
% -----------------------------------------------------------


% Monte Carlo simulation
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- Monte Carlo simulation --- ');
disp(' ');

% initial time of analysis
t0 = 7;

% final time of analysis
% (analysis starting at january 2016)
t1 = 365;

% value of time steps
dt = 1;

% ODE solver optional parameters
%opt = odeset('RelTol',1.0e-6,'AbsTol',1.0e-9);

% interval of analysis
tspan = t0:dt:t1;

% number of time steps
Ndt = length(tspan);

% physical parameters vector
phys_param = [N bH aH yH bV aV dV];

% preallocate memory for system QoIs
MC_SH = zeros(Ns,Ndt);
MC_EH = zeros(Ns,Ndt);
MC_IH = zeros(Ns,Ndt);
MC_RH = zeros(Ns,Ndt);
MC_SV = zeros(Ns,Ndt);
MC_EV = zeros(Ns,Ndt);
MC_IV = zeros(Ns,Ndt);
MC_C  = zeros(Ns,Ndt);

% preallocate memory for MC conv metric
MC_conv = zeros(1,Ns);

for imc=1:Ns
    
    if mod(imc,sqrt(Ns)) == 0
        disp('')
        disp(imc)
    end

    % initial conditions vector
    IH0 = IH0_samples(imc,1);
    SH0 = N-IH0-EH0-RH0;
     IC = [SH0 EH0 IH0 RH0 SV0 EV0 IV0 C0];
    %IC = [SH0(imc) EH0(imc) IH0(imc) RH0(imc) SV0 EV0 IV0 C0];

    % ODE solver Runge-Kutta45
    [time,y] = ode45(@(t,y)rhs_SEIR_SEI(t,y,phys_param),tspan,IC);

    % time series of susceptible humans (number of ppl)
    MC_SH(imc,:) = y(:,1);

    % time series of incubating humans (number of ppl)
    MC_EH(imc,:) = y(:,2);

    % time series of infectious humans (number of ppl)
    MC_IH(imc,:) = y(:,3);

    % time series of recovered humans (number of ppl)
    MC_RH(imc,:) = y(:,4);

    % time series of proportion of susceptible vectors (number of vectors)
    MC_SV(imc,:) = y(:,5);

    % time series of proportion of incubating vectors (number of vectors)
    MC_EV(imc,:) = y(:,6);

    % time series of proportion of infectious vectors (number of vectors)
    MC_IV(imc,:) = y(:,7);

    % time series of cumulative number of infected humans (number of ppl)
    MC_C(imc,:) = y(:,8);
    
    % compute MC convergence metric
    MC_conv(1,imc) = trapz(time,MC_SH(imc,:).^2 + MC_EH(imc,:).^2 + ...
                                MC_IH(imc,:).^2 + MC_RH(imc,:).^2 + ...
                                MC_SV(imc,:).^2 + MC_EV(imc,:).^2 + ...
                                MC_IV(imc,:).^2 +  MC_C(imc,:).^2);
    %MC_conv(1,imc) = trapz(time,MC_IV(imc,:).^2);

end

% compute MC convergence metric
MC_conv = randvar_mc_conv(MC_conv,Ns);

% number of weeks
Nweeks = 52;

% discrete number of new cases
MC_NewCases = zeros(Ns,Nweeks);

MC_NewCases(:,1) = MC_C(:,1);
for i = 2:Nweeks
    MC_NewCases(:,i) = MC_C(:,(7*i-7)/dt+1) - MC_C(:,(7*(i-1)-7)/dt+1);
end

toc
% -----------------------------------------------------------




% compute the statistics
% -----------------------------------------------------------
tic

disp(' ')
disp(' --- compute statistics --- ');
disp(' ');
disp('    ... ');
disp(' ');

% sample average of random processes
SH_smp_avg = mean(MC_SH);
EH_smp_avg = mean(MC_EH);
IH_smp_avg = mean(MC_IH);
RH_smp_avg = mean(MC_RH);
SV_smp_avg = mean(MC_SV);
EV_smp_avg = mean(MC_EV);
IV_smp_avg = mean(MC_IV);
 C_smp_avg = mean(MC_C);
 NewCases_smp_avg = mean(MC_NewCases);

% time average of random processes
SH_time_avg = mean(MC_SH,2);
EH_time_avg = mean(MC_EH,2);
IH_time_avg = mean(MC_IH,2);
RH_time_avg = mean(MC_RH,2);
SV_time_avg = mean(MC_SV,2);
EV_time_avg = mean(MC_EV,2);
IV_time_avg = mean(MC_IV,2);
 C_time_avg = mean(MC_C,2);
NewCases_time_avg = mean(MC_NewCases,2);

% standard deviation of random processes
SH_std = std(MC_SH);
EH_std = std(MC_EH);
IH_std = std(MC_IH);
RH_std = std(MC_RH);
SV_std = std(MC_SV);
EV_std = std(MC_EV);
IV_std = std(MC_IV);
 C_std = std(MC_C);
NewCases_std = std(MC_NewCases);

% skewness of random processes
SH_skew = skewness(MC_SH);
EH_skew = skewness(MC_EH);
IH_skew = skewness(MC_IH);
RH_skew = skewness(MC_RH);
SV_skew = skewness(MC_SV);
EV_skew = skewness(MC_EV);
IV_skew = skewness(MC_IV);
 C_skew = skewness(MC_C);
NewCases_skew = skewness(MC_NewCases);

% kurtosis of the random processes
SH_kurt = kurtosis(MC_SH)-3;
EH_kurt = kurtosis(MC_EH)-3;
IH_kurt = kurtosis(MC_IH)-3;
RH_kurt = kurtosis(MC_RH)-3;
SV_kurt = kurtosis(MC_SV)-3;
EV_kurt = kurtosis(MC_EV)-3;
IV_kurt = kurtosis(MC_IV)-3;
 C_kurt = kurtosis(MC_C)-3;
NewCases_kurt = kurtosis(MC_NewCases)-3;

% confidence probability (percentual)
Pc = 95;

% upper percentil
r_plus = 0.5*(100 + Pc);

% lower percentil
r_minus = 0.5*(100 - Pc);

% confidence bands upper bounds
SH_upp = prctile(MC_SH,r_plus);
EH_upp = prctile(MC_EH,r_plus);
IH_upp = prctile(MC_IH,r_plus);
RH_upp = prctile(MC_RH,r_plus);
SV_upp = prctile(MC_SV,r_plus);
EV_upp = prctile(MC_EV,r_plus);
IV_upp = prctile(MC_IV,r_plus);
 C_upp = prctile( MC_C,r_plus);
NewCases_upp = prctile(MC_NewCases,r_plus);

% confidence bands lower bounds
SH_low = prctile(MC_SH,r_minus);
EH_low = prctile(MC_EH,r_minus);
IH_low = prctile(MC_IH,r_minus);
RH_low = prctile(MC_RH,r_minus);
SV_low = prctile(MC_SV,r_minus);
EV_low = prctile(MC_EV,r_minus);
IV_low = prctile(MC_IV,r_minus);
 C_low = prctile( MC_C,r_minus);
NewCases_low = prctile(MC_NewCases,r_minus);   
 
% number of bins
Nbins = round(sqrt(Ns));

% number of KSD points
Nksd = 100;

% time-average histograms estimation
% [SH_bins,SH_freq] = randvar_pdf(randvar_normalize(SH_time_avg),Nbins);
% [EH_bins,EH_freq] = randvar_pdf(randvar_normalize(EH_time_avg),Nbins);
% [IH_bins,IH_freq] = randvar_pdf(randvar_normalize(IH_time_avg),Nbins);
% [RH_bins,RH_freq] = randvar_pdf(randvar_normalize(RH_time_avg),Nbins);
% [SV_bins,SV_freq] = randvar_pdf(randvar_normalize(SV_time_avg),Nbins);
% [EV_bins,EV_freq] = randvar_pdf(randvar_normalize(EV_time_avg),Nbins);
% [IV_bins,IV_freq] = randvar_pdf(randvar_normalize(IV_time_avg),Nbins);
% [ C_bins, C_freq] = randvar_pdf(randvar_normalize( C_time_avg),Nbins);
% [ NewCases_bins, NewCases_freq] = randvar_pdf(randvar_normalize(NewCases_time_avg),Nbins);

% time-average PDFs curves estimation
% [SH_ksd,SH_supp] = randvar_ksd(randvar_normalize(SH_time_avg),Nksd);
% [EH_ksd,EH_supp] = randvar_ksd(randvar_normalize(EH_time_avg),Nksd);
% [IH_ksd,IH_supp] = randvar_ksd(randvar_normalize(IH_time_avg),Nksd);
% [RH_ksd,RH_supp] = randvar_ksd(randvar_normalize(RH_time_avg),Nksd);
% [SV_ksd,SV_supp] = randvar_ksd(randvar_normalize(SV_time_avg),Nksd);
% [EV_ksd,EV_supp] = randvar_ksd(randvar_normalize(EV_time_avg),Nksd);
% [IV_ksd,IV_supp] = randvar_ksd(randvar_normalize(IV_time_avg),Nksd);
% [ C_ksd, C_supp] = randvar_ksd(randvar_normalize( C_time_avg),Nksd);
% [ NewCases_ksd, NewCases_supp] = randvar_pdf(randvar_normalize(NewCases_time_avg),Nksd);

% random processes marginal histograms estimation
%[SH_bins_t,SH_freq_t] = randvar_pdf(randvar_normalize(MC_SH),Nbins);
% [EH_bins_t,EH_freq_t] = randvar_pdf(randvar_normalize(MC_EH),Nbins);
% [IH_bins_t,IH_freq_t] = randvar_pdf(randvar_normalize(MC_IH),Nbins);
% [RH_bins_t,RH_freq_t] = randvar_pdf(randvar_normalize(MC_RH),Nbins);
% [SV_bins_t,SV_freq_t] = randvar_pdf(randvar_normalize(MC_SV),Nbins);
% [EV_bins_t,EV_freq_t] = randvar_pdf(randvar_normalize(MC_EV),Nbins);
% [IV_bins_t,IV_freq_t] = randvar_pdf(randvar_normalize(MC_IV),Nbins);
% [ C_bins_t, C_freq_t] = randvar_pdf(randvar_normalize( MC_C),Nbins);
% [ NewCases_bins_t, NewCases_freq_t] = randvar_pdf(randvar_normalize(MC_NewCases),Nbins);

% random processes marginal PDFs curves estimation
% [SH_ksd_t,SH_supp_t] = randvar_pdf(randvar_normalize(MC_SH),Nksd);
% [EH_ksd_t,EH_supp_t] = randvar_pdf(randvar_normalize(MC_EH),Nksd);
% [IH_ksd_t,IH_supp_t] = randvar_pdf(randvar_normalize(MC_IH),Nksd);
% [RH_ksd_t,RH_supp_t] = randvar_pdf(randvar_normalize(MC_RH),Nksd);
% [SV_ksd_t,SV_supp_t] = randvar_pdf(randvar_normalize(MC_SV),Nksd);
% [EV_ksd_t,EV_supp_t] = randvar_pdf(randvar_normalize(MC_EV),Nksd);
% [IV_ksd_t,IV_supp_t] = randvar_pdf(randvar_normalize(MC_IV),Nksd);
% [ C_ksd_t, C_supp_t] = randvar_pdf(randvar_normalize( MC_C),Nksd);
% [ NewCases_ksd_t, NewCases_supp_t] = randvar_pdf(randvar_normalize(MC_NewCases),Nksd);

% random variables statistics
SH_time_avg_mean =     mean(SH_time_avg)
SH_time_avg_std  =      std(SH_time_avg)
SH_time_avg_skew = skewness(SH_time_avg)
SH_time_avg_kurt = kurtosis(SH_time_avg)-3

EH_time_avg_mean =     mean(EH_time_avg)
EH_time_avg_std  =      std(EH_time_avg)
EH_time_avg_skew = skewness(EH_time_avg)
EH_time_avg_kurt = kurtosis(EH_time_avg)-3

IH_time_avg_mean =     mean(IH_time_avg)
IH_time_avg_std  =      std(IH_time_avg)
IH_time_avg_skew = skewness(IH_time_avg)
IH_time_avg_kurt = kurtosis(IH_time_avg)-3

RH_time_avg_mean =     mean(RH_time_avg)
RH_time_avg_std  =      std(RH_time_avg)
RH_time_avg_skew = skewness(RH_time_avg)
RH_time_avg_kurt = kurtosis(RH_time_avg)-3

SV_time_avg_mean =     mean(SV_time_avg)
SV_time_avg_std  =      std(SV_time_avg)
SV_time_avg_skew = skewness(SV_time_avg)
SV_time_avg_kurt = kurtosis(SV_time_avg)-3

EV_time_avg_mean =     mean(EV_time_avg)
EV_time_avg_std  =      std(EV_time_avg)
EV_time_avg_skew = skewness(EV_time_avg)
EV_time_avg_kurt = kurtosis(EV_time_avg)-3

IV_time_avg_mean =     mean(IV_time_avg)
IV_time_avg_std  =      std(IV_time_avg)
IV_time_avg_skew = skewness(IV_time_avg)
IV_time_avg_kurt = kurtosis(IV_time_avg)-3

C_time_avg_mean =     mean(C_time_avg)
C_time_avg_std  =      std(C_time_avg)
C_time_avg_skew = skewness(C_time_avg)
C_time_avg_kurt = kurtosis(C_time_avg)-3

NewCases_time_avg_mean =     mean(NewCases_time_avg)
NewCases_time_avg_std  =      std(NewCases_time_avg)
NewCases_time_avg_skew = skewness(NewCases_time_avg)
NewCases_time_avg_kurt = kurtosis(NewCases_time_avg)-3

toc
% -----------------------------------------------------------



% save simulation data into a file
% -----------------------------------------------------------
tic

disp(' ')
disp(' --- Save Workspace --- ');
disp(' ');
disp('    ... ');
disp(' ');


save([num2str(case_name),'.mat']);

toc
% -----------------------------------------------------------


% -----------------------------------------------------------
post_ZikaVD_SEIR_SEI_MonteCarlo
% -----------------------------------------------------------


    
