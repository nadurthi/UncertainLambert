function [Ehat]=keplermod(dt, r0, a, sig0, Ehat0, tol,MU)
% Code to solve the Modified Kepler's Equation
% M (=sqrt(mue/a^3)*dt) - (Eh-(1-r0/a)*sin(E)+sig0*(cos(Eh)-1)/sqrt(a))=0 
% M - Differential Mean anomaly (M-M0)
% sig0 - r0'*v0/sqrt(mue)
% a - semi major axis
% Ehat0 - Eccentric anomaly (initial guess)
% Ehat - Differential Eccentric anomaly (E-E0)
% tol - tolerance for Newton's Root Solver

mue = MU;  % Gravitational parameter of earth

nitr = 100; Ehat = Ehat0;
    M = sqrt(mue/a^3)*dt;

for k = 1:nitr
    Mhat = Ehat - (1-r0/a)*sin(Ehat)-sig0*(cos(Ehat)-1)/sqrt(a);
    rhat = a*(1+(r0/a - 1)*cos(Ehat) + sig0*sin(Ehat)/sqrt(a));
    fx = M - Mhat;
    dfx = -(rhat/a);
    
    dE = -fx/dfx;
    
    if (abs(dE)>tol)
         k; 
         Ehat = Ehat + max(-1,min(dE,1));   % Differential Correction of Newton's Method
    elseif (abs(dE)<tol)
    %    disp(['Solution Converged in ',num2str(k),' iterations']);
        break;
    elseif (abs(dE)>tol)&&(k==nitr)
        disp(['WARNING: Solver Could Not Converge: ',num2str(k),' iterations']);
        break;
    end
end
    