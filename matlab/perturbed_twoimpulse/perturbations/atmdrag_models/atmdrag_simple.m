function adrag = atmdrag_simple(t,x,constants,Satmodel,config)

rsat=x(1:3);
vsat=x(4:6);
% Atmospheric Drag Calculation
%%
Bal = Satmodel.Cd*Satmodel.area/Satmodel.msat; % Ballistic coefficient
R = norm(rsat);

omega_earth = 2*pi/(24*60*60); % Earth's rotation speed in radians per second
vatmos = [-omega_earth*rsat(2); omega_earth*rsat(1);0];   % km/sec
vrel = vsat-vatmos;  % relative velocity of satellite in km/sec
mvrel = norm(vrel);

rho = atmData_simple(R-constants.Re); %kg/m^3
rho = rho/(1e-3)^3; %kg/km^3

adrag =  -1/2*(Bal*rho*mvrel^2)*(vrel/mvrel);
% return;
end