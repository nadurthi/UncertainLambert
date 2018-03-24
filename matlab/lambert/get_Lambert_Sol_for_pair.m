
function [SimSol_Nmax,SimSol_Asol,SimSol_transferorbits]=get_Lambert_Sol_for_pair(XA,XB,t_des,deptime_normalized,constants,config)
XA=XA(:);
XB=XB(:);


r0A = XA(1:3);
v0A = XA(4:6);

r0B = XB(1:3);
v0B = XB(4:6);

r1=norm(r0A);
r2=norm(r0B);

SimSol_Nmax=cell(1,1);
SimSol_Asol=[];
SimSol_transferorbits=[];

theta      = acos(dot(r0A,r0B)/(r1*r2));

% for a given r1,r2,theta ... get lambert problem parameters
lambertparams=GetLambert_parameters(r1,r2,theta,constants,config);

% first find Nmax for the given t_des
Nmax_Output=Get_Nmax(t_des,lambertparams,constants,config);
SimSol_Nmax{1,1}=theta;
SimSol_Nmax{1,2}=Nmax_Output;


% then find 2*Nmax+1 a solutions for the given Nmax
if Nmax_Output.parabolic ==0 && isfinite(Nmax_Output.Nmax) % Parabolic check
    
    %get the 2*Nmax+1 number of a, corresponding branch
    Asol=Get_a(t_des,Nmax_Output,lambertparams,constants,config);
    SimSol_Asol=Asol;
    
    % Solve for the trajectories
    SimSol_transferorbits=Get_all_transferorbit_properties(Asol,r0A,r0B,v0A,v0B,deptime_normalized,t_des,lambertparams,constants,config);
    
end



theta =2*pi-theta; %is also a possibility, next version considers


% for a given r1,r2,theta ... get all the parameters for the
% lambert problem
lambertparams=GetLambert_parameters(r1,r2,theta,constants,config);

% first find Nmax for the given t_des
Nmax_Output=Get_Nmax(t_des,lambertparams,constants,config);
SimSol_Nmax{2,1}=theta;
SimSol_Nmax{2,2}=Nmax_Output;


% then find 2*Nmax+1 a solutions for the given Nmax
if Nmax_Output.parabolic ==0 && isfinite(Nmax_Output.Nmax) % Parabolic check
    
    %get the 2*Nmax+1 number of a, corresponding branch
    Asol=Get_a(t_des,Nmax_Output,lambertparams,constants,config);
    SimSol_Asol=vertcat(SimSol_Asol,Asol);
    
    % Solve for the trajectories
    SimSol_transferorbits=vertcat(SimSol_transferorbits,Get_all_transferorbit_properties(Asol,r0A,r0B,v0A,v0B,deptime_normalized,t_des,lambertparams,constants,config));
    
end


end