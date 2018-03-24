function T=gettimeperiod(r,v,MU)
h=cross(r,v);
magh=norm(h);
a=MU*norm(r)/(2*MU-norm(v)^2);
T=2*pi/sqrt(MU) * a^(3/2);


