% -----------------------------------------------------------------
% ARBO - Arbovirus Modeling and Uncertainty Quantification Toolbox
% -----------------------------------------------------------------

function F = TruncExpParam(lambda,param)

lamb0 = lambda(1);
lamb1 = lambda(2);

xmin = param(1);
xmax = param(2);
mu   = param(3);

F(1) = exp(lamb0-lamb1*xmin)-exp(lamb0-lamb1*xmax)-lamb1;
F(2) = exp(lamb0-lamb1*xmin)*(1.0+lamb1*xmin) - ...
       exp(lamb0-lamb1*xmax)*(1.0+lamb1*xmax) - mu*lamb1^2;
end
