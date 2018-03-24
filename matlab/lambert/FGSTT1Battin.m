function  [T,Y] = FGSTT1Battin(x0,tvec,MU)
% F and G series solution - analytic solution to the 2 Body Problem
% Analytic computation of the state transition tensors
% First order state transition matrix using Schaub and Junkins' Expression
% of Battins' Form 


mue =MU;  % Gravitational parameter of earth
tol = 1e-15; N = 6;
R0 = x0(1:3,:); V0 = x0(4:6,:);
r0 = sqrt(R0'*R0); v0 = sqrt(V0'*V0);
a = 1/(2/r0 - v0^2/mue);
sig0 = R0'*V0/sqrt(mue);
I6 = eye(6); Y(1,1:N+N^2+N^3)=0;
Y(1,1:N)=x0'; T = tvec;
Y(1,N+1:N+N^2) = [I6([1:N^2])]';
Ehat = 0; I3 = eye(3);

for k = 2:length(tvec)
    dt = tvec(k)-tvec(1);
    Ehat0 = Ehat;
    [Ehat]=keplermod(dt, r0, a, sig0, Ehat0, tol,MU);
    Ehatt(k,1)=Ehat;
    
    r = a + (r0 - a)*cos(Ehat) + sqrt(a)*sig0*sin(Ehat);
    % F and G functions
    F(k,1) = 1-a*(1-cos(Ehat))/r0;
    G(k,1) = dt + sqrt(a^3/mue)*(sin(Ehat)-Ehat);
    % Ft and Gt functions
    F(k,2) = -sqrt(mue*a)*sin(Ehat)/(r*r0);
    G(k,2) = 1+a*(cos(Ehat)-1)/r;
    
    Rt = F(k,1)*R0 + G(k,1)*V0;
    Vt = F(k,2)*R0 + G(k,2)*V0;
    Y(k,1:3) = Rt'; 
    Y(k,4:6) = Vt';
    
    %%% State Transition Tensors
    %% First Order
    dr0dR0 = R0/r0; dr0dV0 = 0*R0; 
    dv0dR0 = 0*V0; dv0dV0 = V0/v0;
    dadR0 = 2*a^2*R0/(r0^3); dadV0 = 2*a^2*V0/(mue);
    dsig0dR0 = V0/sqrt(mue); dsig0dV0 = R0/sqrt(mue);
    Ea = (sig0*(1-cos(Ehat))/(2*r*sqrt(a))-1.5*sqrt(mue/a^3)*dt/r +...
        r0*sin(Ehat)/(a*r));
    Er0 = -sin(Ehat)/r; 
    Esig0 = -sqrt(a)*(1-cos(Ehat))/r;
    dEhatdR0 = Ea*dadR0 + Er0*dr0dR0 + Esig0*dsig0dR0;
    dEhatdV0 = Ea*dadV0 + Er0*dr0dV0 + Esig0*dsig0dV0;
    rhoa = 1 - cos(Ehat) + sig0*sin(Ehat)/(2*sqrt(a));
    rhoE = sig0*sqrt(a)*cos(Ehat) - (r0-a)*sin(Ehat);
    rhosig0 = sqrt(a)*sin(Ehat);
    rhor0 = cos(Ehat);
    drdR0 = rhoa*dadR0+rhoE*dEhatdR0+rhor0*dr0dR0+rhosig0*dsig0dR0;
    drdV0 = rhoa*dadV0+rhoE*dEhatdV0+rhor0*dr0dV0+rhosig0*dsig0dV0;
    %%%%% F and G general sensitivity coefficients
    sinE = sin(Ehat); cosE=cos(Ehat);
    Fa=(cosE-1)/r0; Fr0=a*(1-cosE)/(r0^2); FE=-a*sinE/r0;
    dFdR0 = Fa*dadR0 + Fr0*dr0dR0 + FE*dEhatdR0;
    dFdV0 = Fa*dadV0 + Fr0*dr0dV0 + FE*dEhatdV0;
    
    Ga=1.5*sqrt(a/mue)*(sinE-Ehat);GE=sqrt(a^3/mue)*(cosE-1);
    dGdR0 = Ga*dadR0 +  GE*dEhatdR0;
    dGdV0 = Ga*dadV0 +  GE*dEhatdV0;
    
    Fta=-sqrt(mue/a)*sinE/(2*r*r0);Ftr=sqrt(mue*a)*sinE/(r0*r^2);
    Ftr0=sqrt(mue*a)*sinE/(r*r0^2);FtE=-sqrt(mue*a)*cosE/(r0*r);
    dFtdR0 = Fta*dadR0 + Ftr*drdR0 + Ftr0*dr0dR0 + FtE*dEhatdR0;
    dFtdV0 = Fta*dadV0 + Ftr*drdV0 + Ftr0*dr0dV0 + FtE*dEhatdV0;
    
    Gta = -(1-cosE)/r; Gtr=a*(1-cosE)/r^2; GtE=-a*sinE/r;
    dGtdR0 = Gta*dadR0 + Gtr*drdR0 + GtE*dEhatdR0;
    dGtdV0 = Gta*dadV0 + Gtr*drdV0 + GtE*dEhatdV0;

%% Battin's Expression
    dV = Vt-V0; dR=Rt-R0;
    C = a*(3*sinE-(2+cosE)*Ehat)*sqrt(a^3/mue)-dt*a*(1-cosE);
    Phi11=r*dV*dV'/mue+(r0*(1-F(k,1))*Rt+C*Vt)*R0'/r0^3+F(k,1)*I3;
    Phi12=r0*(1-F(k,1))*(dR*V0'-dV*R0')/mue+C*Vt*V0'/mue+G(k,1)*I3;
    Phi21=-dV*R0'/r0^2-Rt*dV'/r^2-mue*C*Rt*R0'/(r*r0)^3;
    Phi21=Phi21+F(k,2)*(I3-Rt*Rt'/r^2+(Rt*Vt'-Vt*Rt')*Rt*dV'/(mue*r));
    Phi22=r0*dV*dV'/mue+(r0*(1-F(k,1))*Rt*R0'-C*Rt*V0')/r^3+G(k,2)*I3;
    
    
%     Phi11 = F(k,1)*I3 + R0*dFdR0' + V0*dGdR0';
%     Phi21 = F(k,2)*I3 + R0*dFtdR0' + V0*dGtdR0';
%     Phi12 = G(k,1)*I3 + R0*dFdV0' + V0*dGdV0';
%     Phi22 = G(k,2)*I3 + R0*dFtdV0' + V0*dGtdV0';
     Phit = [Phi11, Phi12; Phi21, Phi22];
    Y(k,N+1:N+N^2) = [Phit([1:N^2])]';
                    
end

