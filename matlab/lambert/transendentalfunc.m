function [f,dfda,dgda]=transendentalfunc(a,N,alpha,beta,constants,mode)

% Transendental Equation
f     = (6*N*pi + 3*(alpha-beta) - (sin(alpha)-sin(beta)))...
    *(sin(alpha-beta) + sin(alpha)-sin(beta)) - 8*(1-cos(alpha-beta));         % Eq. 5.44

xi    = alpha - beta;
eta   = sin(alpha) - sin(beta);
        
dfda  = ((6*N*pi + 3*xi - eta)*(cos(xi) + cos(alpha)) + ...
            (3-cos(alpha))*(sin(xi) + eta) - 8*sin(xi))*(-1/a*tan(alpha/2)) + ...
            ((6*N*pi + 3*xi - eta)*(-cos(xi) - cos(alpha)) +...
            (-3 - cos(beta))*(sin(xi) + eta) + 8*sin(xi))*(-1/a*tan(beta/2));   % Eq. 5.45

        
% Analytical derivative (dt_da)
dgda = 0.5*sqrt(a/constants.muCan)/(sin(alpha-beta)+(sin(alpha)-sin(beta)))*f;