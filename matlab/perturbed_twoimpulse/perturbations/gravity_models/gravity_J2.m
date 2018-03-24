function ad_J2=gravity_J2(t,x,constants,Satmodel,config)
% t is the Ephemeris Time
% x is eci frame ... typically the J200 frame
% All parameters can be sent in using Satmodel

J2=0.001082616; 

rsat=x(1:3);
r = norm(rsat);
C11=rsat(1)/r;C12=rsat(2)/r;C13=rsat(3)/r;

adx=-1.5*J2*(constants.mu/r^2)*(constants.Re/r)^2*(1-5*C13^2)*C11;
ady=-1.5*J2*(constants.mu/r^2)*(constants.Re/r)^2*(1-5*C13^2)*C12;
adz=-1.5*J2*(constants.mu/r^2)*(constants.Re/r)^2*(3-5*C13^2)*C13;

ad_J2 = [adx;ady;adz];

