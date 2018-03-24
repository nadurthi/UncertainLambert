function transorbit=Get_transferorbit(a,N,r1,r2,branch,t_depart_normalized,t_des,lambertparams,constants,config)
% construct and return the transorbit structure, that contains
%   - the trajectory from r1 to r2 at 100 equi-spaced points between 0
%       and t_des
%   - The altitude at these points
%   - The velocities at r1 and r2, and their corresponding change in
%       velocities
%   - Also all the other parameters such as a,s,c,theta,branch etc that all
%       correspond to this particular transfer orbit/trajectory


R1=norm(r1);
R2=norm(r2);


s=lambertparams.s;
c=lambertparams.c;
theta=lambertparams.theta;

% first find if it is lower or upper for given a
% [alphaL,betaL,~,~,~,~,~,tFL]=Get_alpha_beta(a,s,c,theta,N,t_des,'lower',constants);
% [alphaU,betaU,~,~,~,~,~,tFU]=Get_alpha_beta(a,s,c,theta,N,t_des,'upper',constants);

if strcmp(branch,'lower')
    [alphaL,betaL,~,~,~,~,~,tFL]=Get_alpha_beta(a,s,c,theta,N,t_des,'lower',constants);
    alpha=alphaL;
    beta=betaL;
else
    [alphaU,betaU,~,~,~,~,~,tFU]=Get_alpha_beta(a,s,c,theta,N,t_des,'upper',constants);
    alpha=alphaU;
    beta=betaU;
end

% if abs(tFL-t_des)<abs(tFU-t_des)
%     alpha=alphaL;
%     beta=betaL;
% else
%     alpha=alphaU;
%     beta=betaU;
% end

u1 = r1/R1;                         % Eq. 5.35
u2 = r2/R2;
uc = (r2 - r1)/c;

A = sqrt(constants.muCan/4/a)*cot(alpha/2);   % Eq. 5.37a
B = sqrt(constants.muCan/4/a)*cot(beta/2);    % Eq. 5.37b

v1 = (B+A)*uc + (B-A)*u1;           % Eq. 5.36
v2 = (B+A)*uc - (B-A)*u2;


Nt=100;
time  = linspace(0,t_des,Nt);  % Time vector
rFnG  = zeros(Nt,3);
vFnG  = zeros(Nt,3);
t0    = 0;
flag_ode=0;
for i = 1:length(time)

    [r,v,Ehat] = FnG(t0,time(i),r1,v1,constants.muCan);

    rFnG(i,:)  = r';
    vFnG(i,:)  = v';
end

if norm(rFnG(end,:)'-r2(:))>1e-2
    keyboard
end

transorbit.traj   = [rFnG vFnG];
transorbit.r1=r1;
transorbit.r2=r2;
transorbit.N=N;
transorbit.v1=v1;
transorbit.v2=v2;
transorbit.Ehat=Ehat;
transorbit.a=a;
transorbit.s=s;
transorbit.c=c;
transorbit.a_m=lambertparams.a_m;
transorbit.theta=theta;
transorbit.t_des=t_des;
transorbit.branch=branch;
transorbit.normalized_transfer_timesteps=t_depart_normalized+time;

if config.useSPICE==true
    transorbit.true_Altitude=get_sat_earth_altitude(transorbit.normalized_transfer_timesteps*constants.TU,transorbit.traj(:,1:3)*constants.normX2trueX,constants);
else
    transorbit.true_Altitude=sqrt(sum((transorbit.traj(:,1:3)*constants.normX2trueX).^2,2))-constants.Re;
end