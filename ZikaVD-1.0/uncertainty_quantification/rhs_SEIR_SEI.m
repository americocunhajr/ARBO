%% Right hand side of the SEIR_SEI ZIKAV system
% -----------------------------------------------------------------
% This function defines the system of ODEs for the SEIR-SEI ZIKAV system
%
%   SHdot = -bH*SH.*IV;                                                   rate of susceptible humans 
%                                                                           (number of ppl/days)
% 
%   EHdot = bH*SH.*IV - aH*EH;                                            rate of incubating humans 
%                                                                           (number of ppl/days)
% 
%   IHdot = aH*EH - yh*IH;                                                rate of infectious humans 
%                                                                           (number of ppl/days)
% 
%   RHdot = yh*IH;                                                        rate of recovered humans 
%                                                                           (number of ppl/days)
% 
% 
%   SVdot = dV - bV*SV.*(IH/N) - dV*SV;                                   rate of the proportion of susceptible vectors 
%                                                                           (proportion of vectors/days)
% 
%   EVdot = bV*SV.*(IH/N) - aV*EV - dV*EV;                                rate of the proportion of incubating vectors 
%                                                                           (proportion of vectors/days)
% 
%   IVdot = aV*EV - dV*IV;                                                rate of the proportion of infectious vectors
%                                                                           (proportion of vectors/days)
%
%   Cdot  = aH*EH                                                         rate of the number of cumulative infected humans
%                                                                           (number of ppl/days)
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
%  programmer: Eber Dantas
%              eber.paiva@uerj.br
%              
%              Michel Tosin
%              michel.tosin@uerj.br
%
%  last update: Mar 27, 2017
% -----------------------------------------------------------------


%% Function
% -----------------------------------------------------------------
function ydot = rhs_SEIR_SEI(t,y,phys_param)


%% Physical parameters
    
% phys_param = [N bH aH yH bV aV dV];

N   = phys_param(1);    % human population size (number of ppl)
bH  = phys_param(2);    % vector-to-human transmission rate (ppl/vector)
aH  = phys_param(3);    % incubation human frequency (days^-1)
yH  = phys_param(4);    % infectious human frequency (days^-1)
bV  = phys_param(5);    % human-to-vector transmission rate (vector/ppl)
aV  = phys_param(6);    % incubation vector frequency (days^-1)
dV  = phys_param(7);    % vector lifespan frequency (days^-1)


%% System of ODEs 

% y = [SH EH IH RH SV EV IV C] is the solution vector

ydot = zeros(size(y));


ydot(1) = - bH*y(1).*y(7);                                              % rate of susceptible humans
                                                                        % (number of ppl/days)

ydot(2) = bH*y(1).*y(7) - aH*y(2);                                      % rate of incubating humans
                                                                        %(number of ppl/days)

ydot(3) = aH*y(2) - yH*y(3);                                            % rate of infectious humans
                                                                        % (number of ppl/days)

ydot(4) = yH*y(3);                                                      % rate of recovered humans 
                                                                        % (number of ppl/days)


ydot(5) = dV - bV.*y(5).*(y(3)/N) - dV*y(5);                            % rate of the proportion of susceptible vectors 
                                                                        % (proportion of vectors/days)

ydot(6) = bV*y(5).*(y(3)/N) - aV*y(6) - dV*y(6);                        % rate of the proportion of incubating vectors 
                                                                        % (proportion of vectors/days)

ydot(7) = aV*y(6) - dV*y(7);                                            % rate of the proportion of infectious vectors
                                                                        % (proportion of vectors/days)

ydot(8) = aH*y(2);                                                      % rate of the number of cumulative infectious humans
                                                                        % (number of ppl/days)

end
% -----------------------------------------------------------------
