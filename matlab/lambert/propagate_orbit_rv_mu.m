function rv=propagate_orbit_rv_mu(rv0,Tvec,MU)
% Tvec(1) has to be t0 where r0,v0 is given as initial condition

r0=rv0(1:3);
v0=rv0(4:6);
t0=Tvec(1);

rv=zeros(length(Tvec),6);
for i = 1:length(Tvec)
    [r, v]      = FnG(t0,Tvec(i),r0,v0,MU);
    rv(i,1:3)=r;
    rv(i,4:6)=v;
%     rFnG        = [ rFnG; r' ];
%     vFnG        = [ vFnG; v' ];
end