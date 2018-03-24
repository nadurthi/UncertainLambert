function multistage_lambert_uncertain(XA,XB,t_des)
%% config
config.useSPICE=false;
config.spicepath=NaN;
addpath(genpath('matlab'))

constants= get_constants(config);

%% Setup Arrival and Departure Orbits

deptime_normalized=0;
[SimSol_Nmax,SimSol_Asol,SimSol_transferorbits]=get_Lambert_Sol_for_pair(XA,XB,t_Des,deptime_normalized,constants,config);

%%
close all
Nrequired=0;
Nsol= size(efmoutput.SimSol_transferorbits{1},1);
for i=1:Nsol
    i
    M=efmoutput.SimSol_transferorbits{1}(i,:);
    transorbit=convert_transferorbit_mat2struct(M);
    if transorbit.N == Nrequired
        break
    end
end
[i,transorbit.N,Nsol]
% uncertaintu

rA=transorbit.r1;
rB=transorbit.r2;
v1=transorbit.v1;
v2=transorbit.v2;

t_des = transorbit.t_des;

vA=transorbit.v0A;
vB=transorbit.v0B;

%% nonlinear optimization
close all

Nmc=1000;
PP0=blkdiag((1e-4)^2*eye(3),(1e-4)^2*eye(3));
Xrvmc=mvnrnd([rA,vA],PP0,Nmc);
Wmc=ones(Nmc,1)/Nmc;

Tvec=[0,t_des];

options=optimset('Display','iter','TolX',1e-12);
delv=fminunc(@(delv)optim_ulp(delv,Xrvmc,Wmc,rB,Tvec),v1-vA,options);
% delv=lsqnonlin(@(delv)optim_ulp_lsqnonlin(delv,Xrvmc,Wmc,rB,Tvec),v1-vA,[],[],options);
% delv = fminimax(@(delv)optim_ulp_lsqnonlin(delv,Xrvmc,Wmc,rB,Tvec),v1-vA,[],[],[],[],[],[],[],options);

Nmctest=10000;
Xrvmctest=mvnrnd([rA,vA],1e-2*PP0,Nmctest);
Wmctest=ones(Nmctest,1)/Nmctest;


Err_star=zeros(1,Nmctest);
Err=zeros(1,Nmctest);

dstar=v1-vA;

for i=1:Nmctest
    i
    Xdelv=propagate_orbit_rv_mu([Xrvmctest(i,1),Xrvmctest(i,2),Xrvmctest(i,3),Xrvmctest(i,4)+delv(1),Xrvmctest(i,5)+delv(2),Xrvmctest(i,6)+delv(3)],Tvec,1);
    Xstar=propagate_orbit_rv_mu([Xrvmctest(i,1),Xrvmctest(i,2),Xrvmctest(i,3),Xrvmctest(i,4)+dstar(1),Xrvmctest(i,5)+dstar(2),Xrvmctest(i,6)+dstar(3)],Tvec,1);
    
    Err_star(i) = norm( rB - Xstar(end,1:3) );
    Err(i) = norm( rB - Xdelv(end,1:3) );
    
end