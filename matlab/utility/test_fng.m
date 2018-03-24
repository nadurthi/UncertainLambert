R0     = norm(r0);          % Magnitude of current position
V0     = norm(v0);          % Magnitude of current velocity
sigma0 = (r0'*v0)/sqrt(MU);     % Defined
A      = 2/R0 - V0^2/MU;        % Reciprocal of 1/a
a      = 1/A;                   % Semi-major axis
M      = sqrt(MU/a^3)*(t - t0); % Mean anomaly (rad)

%%
options=optimoptions('fsolve','MaxIter',1e10,'TolFun',1e-12,'Display','off','SpecifyObjectiveGradient',false);
[Ehatfsolve,err]=fsolve(@(x)(M - (x - (1 - R0/a)*sin(x) + sigma0/sqrt(a)*(1 - cos(x)))),12,options);


%%
tol   = 1e-12;       % Tolerance
itr   = 0;          % Initial iteration number
MaxIt = 200;        % Maximum iteration number
Ehat  = M;          % Initial guess for eccentric anomaly (rad)
% if t ~= t0
%     Ehat = Ehat_prev;
% end
% Ehat=0
dEhat = 1;          % Initial eccentric anomaly error (rad)
while abs(dEhat) > tol
    
    err   = M - (Ehat - (1 - R0/a)*sin(Ehat) + sigma0/sqrt(a)*(1 - cos(Ehat)));
    derr  = - 1 + (1 - R0/a)*cos(Ehat) - sigma0/sqrt(a)*sin(Ehat);
    dEhat = max(-1,min(1,err/derr));
    [err,derr,dEhat]
    Ehat  = Ehat -dEhat ;
    
    itr   = itr + 1;
    if itr > MaxIt
        
        disp('hitting max iter for FnG, Switch to alternative method')
      
        
        break
    end
end

Ehat


%%


FEhat   = @(Ehat,M,sigma0,a,R0)M - (Ehat - (1 - R0/a)*sin(Ehat) + sigma0/sqrt(a)*(1 - cos(Ehat)));

Ehat=linspace(0,15*pi,1000);
plot(Ehat,FEhat(Ehat,M,sigma0,a,R0))
grid
