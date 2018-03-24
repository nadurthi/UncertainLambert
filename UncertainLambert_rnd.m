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
Dep_tle_str={'MOLNIYA 2-9',...             
'1 07276U 74026A   18042.27057948  .00000121  00000-0  13303-2 0  9994',...
'2 07276  62.8899 129.6240 6827751 288.5980  12.4440  2.45097420210456'};
% 
% Dep_tle_str={'AGGIESAT 4',...
%     '1 41313U 98067HP  17324.35291052  .00069511  00000-0  26026-3 0  9997',...
%     '2 41313  51.6366 274.9652 0004825  71.5934 288.5595 15.87232926103601'};

Arr_tle_str={'AGGIESAT 4',...
    '1 41313U 98067HP  17324.35291052  .00069511  00000-0  26026-3 0  9997',...
    '2 41313  51.6366 274.9652 0004825  71.5934 288.5595 15.87232926103601'};

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
dept0=DateManager(config).set_date_string('2018-1-21 12:00:00').get_ET();

% Set the Exact end date-time for transfer to begin
deptf=DateManager(config).set_date_string('2017-1-22 1:00:00').get_ET();


% Set the minumum and maximum Transfer time interval
TFt0=0;
TFtf=4*60*60;

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
Nrequired=1;
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
close all
rA=transorbit.r1;
rB=transorbit.r2;
v1=transorbit.v1;
v2=transorbit.v2;

t_des = transorbit.t_des;

vA=transorbit.v0A;
vB=transorbit.v0B;
Nmc=1000;
[xs1,ys1,zs1]= sphere(25);
xs1=xs1*1e-5+rA(1);
ys1=ys1*1e-5+rA(2);
zs1=zs1*1e-5+rA(3);

% [xs2,ys2,zs2]= sphere(50);
% xs2=xs2*1e-6+rA(1);
% ys2=ys2*1e-6+rA(2);
% zs2=zs2*1e-6+rA(3);
%
% [xs3,ys3,zs3]= sphere(50);
% xs3=xs3*1e-7+rA(1);
% ys3=ys3*1e-7+rA(2);
% zs3=zs3*1e-7+rA(3);

Nmc=1000;
Xr1=mvnrnd(rA,1e-6*eye(3),Nmc);

Xtrajmc1=cell(size(xs1,1),size(xs1,2));
% Xtrajmc2=cell(size(xs1,1),size(xs1,2));
% Xtrajmc3=cell(size(xs1,1),size(xs1,2));

Xtrajmc1x_tdes=zeros(size(xs1,1),size(xs1,2));
Xtrajmc1y_tdes=zeros(size(xs1,1),size(xs1,2));
Xtrajmc1z_tdes=zeros(size(xs1,1),size(xs1,2));
%
% Xtrajmc2x_tdes=zeros(size(xs1,1),size(xs1,2));
% Xtrajmc2y_tdes=zeros(size(xs1,1),size(xs1,2));
% Xtrajmc2z_tdes=zeros(size(xs1,1),size(xs1,2));
%
% Xtrajmc3x_tdes=zeros(size(xs1,1),size(xs1,2));
% Xtrajmc3y_tdes=zeros(size(xs1,1),size(xs1,2));
% Xtrajmc3z_tdes=zeros(size(xs1,1),size(xs1,2));


Tvec=[linspace(0,0.0005,50),linspace(0.0005,t_des-2,100),linspace(t_des-2,t_des,100)];
truerv = propagate_orbit_rv_mu([rA,v1],Tvec,1) ;


k=1;
for i=1:size(xs1,1)
    for j=1:size(xs1,2)
        Xtrajmc1{i,j}=propagate_orbit_rv_mu([xs1(i,j),ys1(i,j),zs1(i,j),v1],Tvec,1) ;
        Xtrajmc1x_tdes(i,j)=Xtrajmc1{i,j}(end,1);
        Xtrajmc1y_tdes(i,j)=Xtrajmc1{i,j}(end,2);
        Xtrajmc1z_tdes(i,j)=Xtrajmc1{i,j}(end,3);
        
    end
end

% k=1;
% for i=1:size(xs2,1)
%     for j=1:size(xs2,2)
%         Xtrajmc2{i,j}=propagate_orbit_rv_mu([xs2(i,j),ys2(i,j),zs2(i,j),v1],Tvec,1) ;
%         Xtrajmc2x_tdes(i,j)=Xtrajmc2{i,j}(end,1);
%         Xtrajmc2y_tdes(i,j)=Xtrajmc2{i,j}(end,2);
%         Xtrajmc2z_tdes(i,j)=Xtrajmc2{i,j}(end,3);
%     end
% end
%
%
% k=1;
% for i=1:size(xs3,1)
%     for j=1:size(xs3,2)
%         Xtrajmc3{i,j}=propagate_orbit_rv_mu([xs3(i,j),ys3(i,j),zs3(i,j),v1],Tvec,1) ;
%         Xtrajmc3x_tdes(i,j)=Xtrajmc3{i,j}(end,1);
%         Xtrajmc3y_tdes(i,j)=Xtrajmc3{i,j}(end,2);
%         Xtrajmc3z_tdes(i,j)=Xtrajmc3{i,j}(end,3);
%     end
% end
%

figure(1)
hold on
surface(xs1*constants.normX2trueX,ys1*constants.normX2trueX, zs1*constants.normX2trueX,'FaceColor','green','EdgeColor','none','FaceAlpha',0.4);
camlight right; lighting phong

figure

tri = delaunay(xs1(:),ys1(:));
plot(xs1(:),ys1(:),'.')
[r,c] = size(tri);
disp(r)
h = trisurf(tri, xs1(:),ys1(:), zs1(:),'FaceColor','green','EdgeColor','none');
% l = light('Position',[-50 -15 29])
camlight right; lighting phong
% shading interp
axis equal

% surface(xs2,ys2, zs2,'FaceColor','green','EdgeColor','none','FaceAlpha',0.4);
% camlight right; lighting phong
%
% surface(xs3,ys3, zs3,'FaceColor','green','EdgeColor','none','FaceAlpha',0.4);
% camlight right; lighting phong

axis equal;
plot3(truerv(1:5,1)*constants.normX2trueX,truerv(1:5,2)*constants.normX2trueX,truerv(1:5,3)*constants.normX2trueX,'k','MarkerSize',6,'linewidth',2)

figure(2)
hold on
surface(Xtrajmc1x_tdes*constants.normX2trueX,Xtrajmc1y_tdes*constants.normX2trueX, Xtrajmc1z_tdes*constants.normX2trueX,'FaceColor','green','EdgeColor','none');
%  lighting phong
% axis equal 

% surface(xs1,ys1, zs1,'FaceColor','green','EdgeColor','none','FaceAlpha',0.4);
% camlight right; lighting phong
%
% axis equal;
plot3(truerv(end-2:end,1)*constants.normX2trueX,truerv(end-2:end,2)*constants.normX2trueX,truerv(end-2:end,3)*constants.normX2trueX,'k','MarkerSize',6,'linewidth',2)

figure(3)
hold on
surface(xs1*constants.normX2trueX,ys1*constants.normX2trueX, zs1*constants.normX2trueX,'FaceColor','green','EdgeColor','none','FaceAlpha',0.4);
camlight right; lighting phong
surface(Xtrajmc1x_tdes*constants.normX2trueX,Xtrajmc1y_tdes*constants.normX2trueX, Xtrajmc1z_tdes*constants.normX2trueX,'FaceColor','green','EdgeColor','none','FaceAlpha',0.4);
camlight right; lighting phong
plot3(truerv(:,1)*constants.normX2trueX,truerv(:,2)*constants.normX2trueX,truerv(:,3)*constants.normX2trueX,'k','MarkerSize',6,'linewidth',2)


figure(4)
hold on
plot(Xtrajmc1x_tdes(:),Xtrajmc1y_tdes(:),'ro')

figure(5)
hold on
plot(Xtrajmc1y_tdes(:),Xtrajmc1z_tdes(:),'ro')

figure(6)
hold on
plot(Xtrajmc1x_tdes(:),Xtrajmc1z_tdes(:),'ro')


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


%%
Xtrajmc1_opt=cell(size(xs1,1),size(xs1,2));

Xtrajmc1x_tdes_opt=zeros(size(xs1,1),size(xs1,2));
Xtrajmc1y_tdes_opt=zeros(size(xs1,1),size(xs1,2));
Xtrajmc1z_tdes_opt=zeros(size(xs1,1),size(xs1,2));

k=1;
for i=1:size(xs1,1)
    for j=1:size(xs1,2)
        Xtrajmc1_opt{i,j}=propagate_orbit_rv_mu([xs1(i,j),ys1(i,j),zs1(i,j),v],Tvec,1) ;
        Xtrajmc1x_tdes_opt(i,j)=Xtrajmc1_opt{i,j}(end,1);
        Xtrajmc1y_tdes_opt(i,j)=Xtrajmc1_opt{i,j}(end,2);
        Xtrajmc1z_tdes_opt(i,j)=Xtrajmc1_opt{i,j}(end,3);
        
    end
end

figure(2)

surface(Xtrajmc1x_tdes_opt*constants.normX2trueX,Xtrajmc1y_tdes_opt*constants.normX2trueX, Xtrajmc1z_tdes_opt*constants.normX2trueX,'FaceColor','red','EdgeColor','none','FaceAlpha',0.4);
camlight right; lighting phong

figure(4)
hold on
plot(Xtrajmc1x_tdes_opt(:),Xtrajmc1y_tdes_opt(:),'b*')

figure(5)
hold on
plot(Xtrajmc1y_tdes_opt(:),Xtrajmc1z_tdes_opt(:),'b*')

figure(6)
hold on
plot(Xtrajmc1x_tdes_opt(:),Xtrajmc1z_tdes_opt(:),'b*')

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
Xls=mvnrnd(rA,(1e-5/3)^2*eye(3),Nls);
wls=ones(1,Nls)/Nls;

% xs1,ys1,zs1
PHIls=cell(1,Nls);
Vls=cell(1,Nls);

parfor i=1:Nls
        disp(['i = ',num2str(i)])
        XA=[Xls(i,:),vA];
        XB=[rB,vB];
        deptime_normalized=0;
        [SimSol_Nmax,SimSol_Asol,SimSol_transferorbits]=get_Lambert_Sol_for_pair(XA,XB,t_des,deptime_normalized,constants,config);
        
        Nsol= size(SimSol_transferorbits,1);
        for k=1:Nsol
            k
            M=SimSol_transferorbits(k,:);
            TOLS=convert_transferorbit_mat2struct(M);
            if TOLS.N == Nrequired
                break
            end
        end
        
        [T,Y] = FGSTT1Battin([Xls(i,:),TOLS.v1]',[0,t_des],1);
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

