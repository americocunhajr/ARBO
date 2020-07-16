%% Right hand side of the SEIR_SEI ZIKAV system for the TRR algorithm
% -----------------------------------------------------------------
% This function defines the system of ODEs of the SEIR-SEI ZIKAV system for the TRR algorithm.
% This function is called by the TRR_FunctionOutput_SEIR_SEI to calculate the system response 
% according to the case_name specified
%
%   case_name      = case configuration. Need to be passed down to the auxiliary function 
%   variable_param = parameters varied by the algorithm
%   fixed_ param   = parameters fixed
%   IC             = vector of initial conditions
%
%
%   SHdot = -bH*SH.*IV;                           rate of susceptible humans (number of ppl/days)
%   EHdot = bH*SH.*IV - aH*EH;                    rate of incubating humans (number of ppl/days)
%   IHdot = aH*EH - yh*IH;                        rate of infectious humans number of ppl/days)
%   RHdot = yh*IH;                                rate of recovered humans (number of ppl/days)
% 
%   SVdot = dV - bV*SV.*(IH/N) - dV*SV;           rate of the proportion of susceptible vectors (proportion of vectors/days)
%   EVdot = bV*SV.*(IH/N) - aV*EV - dV*EV;        rate of the proportion of incubating vectors (proportion of vectors/days)
%   IVdot = aV*EV - dV*IV;                        rate of the proportion of infectious vectors(proportion of vectors/days)
%   Cdot  = aH*EH                                 rate of the number of cumulative infected humans (number of ppl/days)
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
%  last update: May 2, 2017
% -----------------------------------------------------------------

%% Function
% -----------------------------------------------------------------
function ydot = TRR_rhs_SEIR_SEI(t,y,case_name,variable_TRRparam,fixed_TRRparam)

%% Case configuration
% -----------------------------------------------------------------

% case_name include the fixed parameters. The variable parameters are specified. 
    % Maintain default parameter vector order!

% Default TRRparam vector: [N bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 IV0 C0]
    % RH0 is omitted to optimize the procedure

if strcmp(case_name,'case_config1');

    %fixed_TRRparam = [N bH aH yH bV aV dV SH0 EH0 SV0 EV0 C0];  
    %inital_guess_TRRparam = [IV0];
    
    N   = fixed_TRRparam(1);  % human population size (number of ppl)
    bH  = fixed_TRRparam(2);  % vector-to-human transmission rate (ppl/vector)    
    aH  = fixed_TRRparam(3);  % incubation human frequency (days^-1)
    yH  = fixed_TRRparam(4);  % infectious human frequency (days^-1)
    bV  = fixed_TRRparam(5);  % human-to-vector transmission rate (vector/ppl)
    aV  = fixed_TRRparam(6);  % incubation vector frequency (days^-1)
    dV  = fixed_TRRparam(7);  % vector lifespan frequency (days^-1)
           
    
elseif strcmp(case_name,'case_config2');
 
    % fixed_TRRparam = [N bH aH yH bV aV dV SH0 EH0 IH0 SV0 C0];  
    % variable_TRRparam = [EV0 IV0];

    N   = fixed_TRRparam(1);  % human population size (number of ppl)
    bH  = fixed_TRRparam(2);  % vector-to-human transmission rate (ppl/vector)    
    aH  = fixed_TRRparam(3);  % incubation human frequency (days^-1)
    yH  = fixed_TRRparam(4);  % infectious human frequency (days^-1)
    bV  = fixed_TRRparam(5);  % human-to-vector transmission rate (vector/ppl)
    aV  = fixed_TRRparam(6);  % incubation vector frequency (days^-1)
    dV  = fixed_TRRparam(7);  % vector lifespan frequency (days^-1)

elseif strcmp(case_name,'case_config3');
 
    % fixed_TRRparam = [N aH yH aV dV SH0 EH0 IH0 SV0 C0];  
    % variable_TRRparam = [bH bV EV0 IV0];

    N   = fixed_TRRparam(1);  % human population size (number of ppl)
    bH  = variable_TRRparam(1);  % vector-to-human transmission rate (ppl/vector)    
    aH  = fixed_TRRparam(2);  % incubation human frequency (days^-1)
    yH  = fixed_TRRparam(3);  % infectious human frequency (days^-1)
    bV  = variable_TRRparam(2);  % human-to-vector transmission rate (vector/ppl)
    aV  = fixed_TRRparam(4);  % incubation vector frequency (days^-1)
    dV  = fixed_TRRparam(5);  % vector lifespan frequency (days^-1)
    
elseif strcmp(case_name,'case_config4');
    
    % fixed_TRRparam = [N C0];
    % variable_TRRparam = [bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 IV0];
    
    N   = fixed_TRRparam(1);     % human population size (number of ppl)
    bH  = variable_TRRparam(1);  % vector-to-human transmission rate (ppl/vector)    
    aH  = variable_TRRparam(2);  % incubation human frequency (days^-1)
    yH  = variable_TRRparam(3);  % infectious human frequency (days^-1)
    bV  = variable_TRRparam(4);  % human-to-vector transmission rate (vector/ppl)
    aV  = variable_TRRparam(5);  % incubation vector frequency (days^-1)
    dV  = variable_TRRparam(6);  % vector lifespan frequency (days^-1)

elseif strcmp(case_name,'case_config5');
    
    % Case with Sums
    % fixed_TRRparam = [N C0];
    % variable_TRRparam = [bH aH yH bV aV dV SH0 EH0 IH0 SV0 EV0 IV0];
    
    N   = fixed_TRRparam(1);     % human population size (number of ppl)
    bH  = variable_TRRparam(1);  % vector-to-human transmission rate (ppl/vector)    
    aH  = variable_TRRparam(2);  % incubation human frequency (days^-1)
    yH  = variable_TRRparam(3);  % infectious human frequency (days^-1)
    bV  = variable_TRRparam(4);  % human-to-vector transmission rate (vector/ppl)
    aV  = variable_TRRparam(5);  % incubation vector frequency (days^-1)
    dV  = variable_TRRparam(6);  % vector lifespan frequency (days^-1)

    
end


%% System of ODEs 
% -----------------------------------------------------------------

% y = [SH EH IH SV EV IV C] is the solution vector
    % RH is omitted to optimize the procedure
    
ydot = zeros(size(y));


ydot(1) = -bH*y(1).*y(6);                           % rate of susceptible humans
                                                        % (number of ppl/days)

ydot(2) = bH*y(1).*y(6) - aH*y(2);                  % rate of incubating humans
                                                        %(number of ppl/days)

ydot(3) = aH*y(2) - yH*y(3);                        % rate of infectious humans
                                                        % (number of ppl/days)

ydot(4) = dV - bV.*y(4).*(y(3)/N) - dV*y(4);        % rate of the proportion of susceptible vectors 
                                                        % (proportion of vectors/days)

ydot(5) = bV*y(4).*(y(3)/N) - aV*y(5) - dV*y(5);    % rate of the proportion of incubating vectors 
                                                        % (proportion of vectors/days)

ydot(6) = aV*y(5) - dV*y(6);                        % rate of the proportion of infectious vectors
                                                        % (proportion of vectors/days)

ydot(7) = aH*y(2);                                  % cumulative number of infectious humans 
                                                        % (number of ppl at t)

end
% -----------------------------------------------------------------
