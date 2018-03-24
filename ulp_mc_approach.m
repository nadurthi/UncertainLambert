%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%% config
config.useSPICE=false;
config.spicepath=NaN;
addpath(genpath('matlab'))

constants= get_constants(config);

%% Setup Arrival and Departure Orbits


% Provide Two Line Elements   %%%%%%%%%%%%%%%%%
% Dep_tle_str={'MOLNIYA 2-9',...             
% '1 07276U 74026A   18042.27057948  .00000121  00000-0  13303-2 0  9994',...
% '2 07276  62.8899 129.6240 6827751 288.5980  12.4440  2.45097420210456'};


Dep_tle_str={'IRIDIUM 10',...         
'1 24839U 97030D   18082.50437020  .00000074  00000-0  19257-4 0  9994',...
'2 24839  86.3961 233.2642 0002034  92.3438 267.7990 14.34212966 86991'};

Arr_tle_str={'IRIDIUM 13',...          
'1 24840U 97030E   18082.47639146  .01265173 -10400-5  11078-2 0  9991',...
'2 24840  86.3947 233.1574 0409371 352.2988   7.2089 15.37286310 87044'};

% Dep_tle_str={'GOES 12 [-]',...             
% '1 26871U 01031A   18059.64693949 -.00000135  00000-0  00000-0 0  9992',...
% '2 26871   6.7859  58.1789 0011524 243.5729 100.4046  0.99139659 60714'};


% Dep_tle_str={'AGGIESAT 4',...
%     '1 41313U 98067HP  17324.35291052  .00069511  00000-0  26026-3 0  9997',...
%     '2 41313  51.6366 274.9652 0004825  71.5934 288.5595 15.87232926103601'};

% Arr_tle_str={'AGGIESAT 4',...
%     '1 41313U 98067HP  17324.35291052  .00069511  00000-0  26026-3 0  9997',...
%     '2 41313  51.6366 274.9652 0004825  71.5934 288.5595 15.87232926103601'};

% Arr_tle_str={'GOES 12 [-]',...             
% '1 26871U 01031A   18059.64693949 -.00000135  00000-0  00000-0 0  9992',...
% '2 26871   6.7859  58.1789 0011524 243.5729 100.4046  0.99139659 60714'}

% 
% Arr_tle_str={'MOLNIYA 2-9',...             
% '1 07276U 74026A   18042.27057948  .00000121  00000-0  13303-2 0  9994',...
% '2 07276  62.8899 129.6240 6827751 288.5980  12.4440  2.45097420210456'};

% Arr_tle_str={'Target: COSMOS 2425 (716)',...
%     '1 29670U 06062A   17324.35552733 -.00000028  00000-0  10000-3 0  9993',...
%     '2 29670  65.3845 299.5924 0021533 338.2090  83.6349  2.13104919 84865'};

% set departure orbit
DepIM=OrbitInputManager(constants);
DepIM.set_TLE(Dep_tle_str,config);
Dep_orb_elem=DepIM.get_orb_elem_vec();
dep_epoch=DepIM.epochDM.get_ET();

% set Arrival orbit
ArrIM=OrbitInputManager(constants);
ArrIM.set_TLE(Arr_tle_str,config);
Arr_orb_elem=ArrIM.get_orb_elem_vec();
arr_epoch=ArrIM.epochDM.get_ET();


depgridlen=1;
TFgridlen=1;

% Set the Exact start date-time for transfer to begin
dept0=DateManager(config).set_date_string('2018-1-19 1:00:00').get_ET();

% Set the Exact end date-time for transfer to begin
deptf=dept0;


% Set the minumum and maximum Transfer time interval
TFt0=0;
TFtf=10*60*60;


% dept0=DateManager(config).set_date_string('2018-1-21 12:00:00').get_ET();
% deptf=DateManager(config).set_date_string('2017-1-22 1:00:00').get_ET();
% TFt0=0;
% TFtf=10*60*60;

disp('done')
% Compute the EFM

disp('Now running EFM')
savefilename='saveddata/Simulation';
efmoutput=EFMcompute_kepler(savefilename,Dep_orb_elem,dep_epoch,Arr_orb_elem,arr_epoch,depgridlen,TFgridlen,dept0,deptf,TFt0,TFtf,constants,config);

%
XA=efmoutput.inputconfig.DepOrbit_normalized;
XB=efmoutput.inputconfig.ArrOrbit_normalized_cell{1}{1};
t_Des=efmoutput.inputconfig.TimeFlights_normalized;
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
    transorbit.N 
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

[xs1,ys1,zs1]=sphere(50);
xs1=xs1*1e-5+rA(1);
ys1=ys1*1e-5+rA(2);
zs1=zs1*1e-5+rA(3);

[vxs1,vys1,vzs1]=sphere(50);

vxs1=vxs1*1e-4+vA(1);
vys1=vys1*1e-4+vA(2);
vzs1=vzs1*1e-4+vA(3);

Nmc=prod(size(xs1));
% Xr1=mvnrnd(rA,1e-6*eye(3),Nmc);

Xtrajmc1=cell(size(xs1,1),size(xs1,2));
Xtrajmc2=cell(size(xs1,1),size(xs1,2));
Xtrajmc3=cell(size(xs1,1),size(xs1,2));

Xtrajmc1x_tdes=zeros(size(xs1,1),size(xs1,2));
Xtrajmc1y_tdes=zeros(size(xs1,1),size(xs1,2));
Xtrajmc1z_tdes=zeros(size(xs1,1),size(xs1,2));
%
Xtrajmc2x_tdes=zeros(size(xs1,1),size(xs1,2));
Xtrajmc2y_tdes=zeros(size(xs1,1),size(xs1,2));
Xtrajmc2z_tdes=zeros(size(xs1,1),size(xs1,2));

Xtrajmc3x_tdes=zeros(size(xs1,1),size(xs1,2));
Xtrajmc3y_tdes=zeros(size(xs1,1),size(xs1,2));
Xtrajmc3z_tdes=zeros(size(xs1,1),size(xs1,2));


Tvec=linspace(0,t_des,5);
% Tvec=linspace(0,t_des,50),linspace(0.0005,t_des-2,100),linspace(t_des-2,t_des,100)];
truerv = propagate_orbit_rv_mu([rA,v1],linspace(0,t_des,100000),1) ;
% truerv = propagate_orbit_rv_mu([rA,v1],linspace(0,t_des,5000),1) ;

TA=gettimeperiod(rA,vA,1);
TB=gettimeperiod(rB,vB,1);
XrvA = propagate_orbit_rv_mu([rA,vA],linspace(0,4*TA,5000),1) ;
XrvB = propagate_orbit_rv_mu([rB,vB],linspace(0,4*TB,500),1) ;



delv=v1-vA;

k=1;
for i=1:size(xs1,1)
    for j=1:size(xs1,2)
        [i,j]

        Xtrajmc1{i}=propagate_orbit_rv_mu([xs1(i,j),ys1(i,j),zs1(i,j),vxs1(i,j)+delv(1),vys1(i,j)+delv(2),vzs1(i,j)+delv(3)],Tvec,1) ;
        Xtrajmc1x_tdes(i,j)=Xtrajmc1{i}(end,1);
        Xtrajmc1y_tdes(i,j)=Xtrajmc1{i}(end,2);
        Xtrajmc1z_tdes(i,j)=Xtrajmc1{i}(end,3);
        
    end
end


figure
hold on
surf(xs1*constants.normX2trueX,ys1*constants.normX2trueX, zs1*constants.normX2trueX,'FaceColor','green','EdgeColor','none','FaceAlpha',0.4);
camlight right; lighting phong
axis equal;
alpha 0.4
plot3(truerv(1:2,1)*constants.normX2trueX,truerv(1:2,2)*constants.normX2trueX,truerv(1:2,3)*constants.normX2trueX,'k','MarkerSize',6,'linewidth',2)


figure
hold on
surf(Xtrajmc1x_tdes*constants.normX2trueX,Xtrajmc1y_tdes*constants.normX2trueX, Xtrajmc1z_tdes*constants.normX2trueX,'FaceColor','green','EdgeColor','none','FaceAlpha',0.4);
camlight right; lighting phong
axis equal 
plot3(truerv(end-105:end,1)*constants.normX2trueX,truerv(end-105:end,2)*constants.normX2trueX,truerv(end-105:end,3)*constants.normX2trueX,'k','MarkerSize',6,'linewidth',2)

figure
plot3(truerv(:,1)*constants.normX2trueX,truerv(:,2)*constants.normX2trueX,truerv(:,3)*constants.normX2trueX,'k','MarkerSize',6,'linewidth',2)
hold on
plot3(XrvA(:,1)*constants.normX2trueX,XrvA(:,2)*constants.normX2trueX,XrvA(:,3)*constants.normX2trueX,'r--','MarkerSize',6,'linewidth',2)
plot3(XrvB(:,1)*constants.normX2trueX,XrvB(:,2)*constants.normX2trueX,XrvB(:,3)*constants.normX2trueX,'b--','MarkerSize',6,'linewidth',2)

ToPo=load('topo.mat','topo','topomap1');
Xsphere.radii=constants.radii;

[x,y,z] = sphere(50);
Xsphere.x=-x;
Xsphere.y=-y;
Xsphere.z=z;
plot_earth_sphere(Xsphere,ToPo)
axis equal
error('end sim')
%% optimize robust Lambert Problem



% xs1,ys1,zs1
PHIs=cell(1,size(xs1,1));
Vs=cell(1,size(xs1,1));
p=1;
parfor i=1:size(xs1,1)
    PHIs{i}=cell(1,size(xs1,2));
    Vs{i}=cell(1,size(xs1,2));
    for j=1:size(xs1,2)
        [i,j]
        XA=[xs1(i,j),ys1(i,j),zs1(i,j),vA];
        XB=[rB,vB];
        deptime_normalized=0;
        [SimSol_Nmax,SimSol_Asol,SimSol_transferorbits]=get_Lambert_Sol_for_pair(XA,XB,t_des,deptime_normalized,constants,config);
        
        Nsol= size(SimSol_transferorbits,1);
        for k=1:Nsol
            k
            M=SimSol_transferorbits(k,:);
            TO=convert_transferorbit_mat2struct(M);
            if TO.N == Nrequired
                break
            end
        end
        
        [T,Y] = FGSTT1Battin([xs1(i,j),ys1(i,j),zs1(i,j),TO.v1]',[0,t_des],1);
        S=reshape(Y(end,1+6:36+6),6,6);
        PHIs{i}{j}=S;
        Vs{i}{j}=TO.v1;
        p=p+1;
    end
end

PHIs2=cell(1,prod(size(xs1)));
Vs2=cell(1,prod(size(xs1)));
p=1;
for i=1:size(xs1,1)
    for j=1:size(xs1,2)
        PHIs2{p}=PHIs{i}{j};
        Vs2{p}=Vs{i}{j};
        p=p+1;
    end
end
%%
[v,ep]=linear_robust_lambert(PHIs2,Vs2);
v=v(:)';

delv_star=v1-vA
delv_rbst=

%%
Xtrajmc1_opt=cell(size(xs1,1),size(xs1,2));
Nmctest=1000;
Xtest0=mvnrnd([rA,vA],1e-5*eye(6),Nmctest);
XtestF_opt=zeros(size(Xtest0));

k=1;
for i=1:Nmctest

        X=propagate_orbit_rv_mu([Xtest0(i,1),Xtest0(i,2),Xtest0(i,3),v],Tvec,1) ;
        XtestF_opt(i,:)=X(end,:);
        

end

%%
Errx=zeros(size(xs1,1),size(xs1,2));
Erry=zeros(size(xs1,1),size(xs1,2));
Errz=zeros(size(xs1,1),size(xs1,2));

for i=1:size(xs1,1)
    for j=1:size(xs1,2)
        
        Errx(i,j)=Xtrajmc1x_tdes_opt(i,j)-Xtrajmc1x_tdes(i,j);
        Erry(i,j)=Xtrajmc1y_tdes_opt(i,j)-Xtrajmc1y_tdes(i,j);
        Errz(i,j)=Xtrajmc1z_tdes_opt(i,j)-Xtrajmc1z_tdes(i,j);
        
    end
end

[min(min(Errx))*constants.normX2trueX,max(max(Errx))*constants.normX2trueX]
[min(min(Erry))*constants.normX2trueX,max(max(Erry))*constants.normX2trueX]
[min(min(Errz))*constants.normX2trueX,max(max(Errz))*constants.normX2trueX]



%% least squares type
Nls=2500;
Xls=mvnrnd([rA,vA],(1e-5)^2*eye(6),Nls);
wls=ones(1,Nls)/Nls;

% xs1,ys1,zs1
PHIls=cell(1,Nls);
Vls=cell(1,Nls);

parfor i=1:Nls
        disp(['i = ',num2str(i)])
        XA=[Xls(i,:),vA];
        [T,Y] = FGSTT1Battin([Xls(i,1:3),v1]',[0,t_des],1);
        S=reshape(Y(end,1+6:36+6),6,6);
        PHIls{i}=S;
        Vls{i}=TOLS.v1;

end

A=0;
B=0;
C=0;

for i=1:Nls
    Si=PHIls{i}(1:3,4:6);
    vi=Vls{i};
    vi=vi(:);
    A = A+ wls(i)*Si'*Si;
    B = B+wls(i)*Si'*Si*vi;
    C = C+wls(i)*vi'*Si'*Si*vi;
    
end

vls = A\B;
vls=vls';
%%
Xtrajmc1_optls=cell(size(xs1,1),size(xs1,2));

Xtrajmc1x_tdes_optls=zeros(size(xs1,1),size(xs1,2));
Xtrajmc1y_tdes_optls=zeros(size(xs1,1),size(xs1,2));
Xtrajmc1z_tdes_optls=zeros(size(xs1,1),size(xs1,2));

k=1;
for i=1:size(xs1,1)
    for j=1:size(xs1,2)
        Xtrajmc1_optls{i,j}=propagate_orbit_rv_mu([xs1(i,j),ys1(i,j),zs1(i,j),vls],Tvec,1) ;
        Xtrajmc1x_tdes_optls(i,j)=Xtrajmc1_optls{i,j}(end,1);
        Xtrajmc1y_tdes_optls(i,j)=Xtrajmc1_optls{i,j}(end,2);
        Xtrajmc1z_tdes_optls(i,j)=Xtrajmc1_optls{i,j}(end,3);
        
    end
end

figure(2)

surface(Xtrajmc1x_tdes_optls*constants.normX2trueX+180,Xtrajmc1y_tdes_optls*constants.normX2trueX, Xtrajmc1z_tdes_optls*constants.normX2trueX,'FaceColor','red','EdgeColor','none','FaceAlpha',0.4);
camlight right; lighting phong


figure(4)
hold on
plot(Xtrajmc1x_tdes_optls(:),Xtrajmc1y_tdes_optls(:),'gs')

figure(5)
hold on
plot(Xtrajmc1y_tdes_optls(:),Xtrajmc1z_tdes_optls(:),'gs')

figure(6)
hold on
plot(Xtrajmc1x_tdes_optls(:),Xtrajmc1z_tdes_optls(:),'gs')

%%
Errxoptls=zeros(size(xs1,1),size(xs1,2));
Erryoptls=zeros(size(xs1,1),size(xs1,2));
Errzoptls=zeros(size(xs1,1),size(xs1,2));

for i=1:size(xs1,1)
    for j=1:size(xs1,2)
        
        Errxoptls(i,j)=Xtrajmc1x_tdes_optls(i,j)-Xtrajmc1x_tdes(i,j);
        Erryoptls(i,j)=Xtrajmc1y_tdes_optls(i,j)-Xtrajmc1y_tdes(i,j);
        Errzoptls(i,j)=Xtrajmc1z_tdes_optls(i,j)-Xtrajmc1z_tdes(i,j);
        
    end
end

[min(min(Errxoptls))*constants.normX2trueX,max(max(Errxoptls))*constants.normX2trueX]
[min(min(Erryoptls))*constants.normX2trueX,max(max(Erryoptls))*constants.normX2trueX]
[min(min(Errzoptls))*constants.normX2trueX,max(max(Errzoptls))*constants.normX2trueX]

%%

dopt=zeros(prod(size(xs1)));
doptls=zeros(prod(size(xs1)));

dopt= sum([Xtrajmc1x_tdes_opt(:)-Xtrajmc1x_tdes(:),Xtrajmc1y_tdes_opt(:)-Xtrajmc1y_tdes(:),Xtrajmc1z_tdes_opt(:)-Xtrajmc1z_tdes(:)].^2,2);
doptls= sum([Xtrajmc1x_tdes_optls(:)-Xtrajmc1x_tdes(:),Xtrajmc1y_tdes_optls(:)-Xtrajmc1y_tdes(:),Xtrajmc1z_tdes_optls(:)-Xtrajmc1z_tdes(:)].^2,2);

[dopt,doptls]

dopt= sqrt(sum([Xtrajmc1x_tdes_opt(:)-rB(1),Xtrajmc1y_tdes_opt(:)-rB(2),Xtrajmc1z_tdes_opt(:)-rB(3)].^2,2))*constants.normX2trueX;
doptls= sqrt(sum([Xtrajmc1x_tdes_optls(:)-rB(1),Xtrajmc1y_tdes_optls(:)-rB(2),Xtrajmc1z_tdes_optls(:)-rB(3)].^2,2))*constants.normX2trueX;
dgrid= sqrt(sum([Xtrajmc1x_tdes(:)-rB(1),Xtrajmc1y_tdes(:)-rB(2),Xtrajmc1z_tdes(:)-rB(3)].^2,2))*constants.normX2trueX;
[dopt,doptls,dgrid]

figure
histogram(dopt,linspace(0,500,200),'Normalization','cdf','FaceColor','r','FaceAlpha',0.4)
hold on
histogram(dgrid,linspace(0,500,200),'Normalization','cdf','FaceColor','b','FaceAlpha',0.4)

% axis([0,5,0,70])

figure
histogram(doptls,linspace(0,500,200),'Normalization','cdf','FaceColor','r','FaceAlpha',0.4)
hold on
histogram(dgrid,linspace(0,500,200),'Normalization','cdf','FaceColor','b','FaceAlpha',0.4)

% hist(doptls,linspace(0,5,100),'r')
% axis([0,5,0,70])

figure
histogram(dgrid,linspace(0,500,200),'Normalization','cdf')
% hist(dgrid,linspace(0,5,100),'r')
% axis([0,5,0,70])

[max(dopt),max(doptls),max(dgrid)]
[mean(dopt),mean(doptls),mean(dgrid)]
[std(dopt),std(doptls),std(dgrid)]

