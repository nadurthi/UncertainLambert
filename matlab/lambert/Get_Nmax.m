function output=Get_Nmax(t_des,lambertparams,constants,mode)
% get the Nmax for the given t_des
% this is done by iterating over N and finding the largest N for which
% t_des > t_minN.
% where t_minN is the minimum transfer time for a given N

% output contains the Nmax and Min_energy_Min_Time_Sols
% Min_energy_Min_Time_Sols: is a matrix that contains for every N<=N_max
%     1. a_m : the a corresponding to minimum energy transfer
%     2. t_m : transfer time corresponding to minimum energy
%     3. a_min: the a corresponding to the smallest trasnfer time for a
%     given N
%     4. t_min : smallest time of transfer for given N
maxiter=150;

% keyboard

theta=lambertparams.theta;
s=lambertparams.s;
c=lambertparams.c;

%% Check if parabolic orbit

if t_des < lambertparams.t_p
    Min_energy_Min_Time_Sols=[NaN,NaN,NaN,NaN,NaN,t_des];%[N,a_m,t_mN,a_min,t_min,t_des]
    output.Nmax=NaN;
    output.Min_energy_Min_Time_Sols=Min_energy_Min_Time_Sols;
    output.parabolic=1;
    
    return
end

%% Check if Nmax=0 is the only possibility
if t_des <= lambertparams.t_m && t_des > lambertparams.t_p
    Min_energy_Min_Time_Sols=[0,lambertparams.a_m,lambertparams.t_m,NaN,NaN,t_des];%[N,a_m,t_mN,a_min,t_min,t_des]
    output.Nmax=0;
    output.Min_energy_Min_Time_Sols=Min_energy_Min_Time_Sols;
    output.parabolic=0;
    
    return
end
%% Determine Nmax
N= 0;  % Starting N value
tF = 0;
Min_energy_Min_Time_Sols=[0,lambertparams.a_m,lambertparams.t_m,NaN,NaN,t_des]; %[N,a_m,t_mN,a_min,t_min,t_des]
P=0;
while tF < t_des
    Nmax = N;
    N    = N + 1;
    f    = 10;
    itr  = 0;
    
    a = 1.001*lambertparams.a_m; % Initial guess for semimajor axis
    [alpha,beta,~,~,a_m,~,t_mN,tF]=Get_alpha_beta(a,s,c,theta,N,t_des,'lower',constants,mode);
    
    if N>=13    
%         keyboard
        itr=1;
        aL=lambertparams.a_m;
        aU= 1.001*aL;
        a=(aL+aU)/2;
        % first get the bounds around the t_min
        cnt=1;
        while true
            
            [~,~,~,~,~,~,~,tL]=Get_alpha_beta(aL,s,c,theta,N,t_des,'lower',constants,mode);
            [~,~,~,~,~,~,~,tU]=Get_alpha_beta(aU,s,c,theta,N,t_des,'lower',constants,mode);
            [~,~,~,~,~,~,~,tmid]=Get_alpha_beta(a,s,c,theta,N,t_des,'lower',constants,mode);
            if tmid<=tU && tmid <= tL
                break
            end
            aL=aU;
            aU= 1.001*aL;
            a=(aL+aU)/2;
            cnt=cnt+1;
        end
        err=1;
        while(err)>1e-10
            [alphamid,betamid,~,~,a_m,~,t_mN,tF]=Get_alpha_beta(a,s,c,theta,N,t_des,'lower',constants,mode);
            [fmid,~,~]=transendentalfunc(a,N,alphamid,betamid,constants,mode);
            
            if fmid>0
                aU=a;
            else
                aL=a;
            end
            a=(aU+aL)/2;

            err=abs(fmid);
            
            itr=itr+1;
            if itr>50
                disp('Using fsolve as backup')
                options=optimoptions('fsolve','MaxIter',1e10,'TolFun',1e-10,'Display','off','SpecifyObjectiveGradient',true);
                [a,err]=fsolve(@(a)solve_for_amin(a,N,alpha,beta,constants,mode),a,options);
                [alpha,beta,~,~,a_m,~,t_mN,tF]=Get_alpha_beta(a,s,c,theta,N,t_des,'lower',constants,mode);
                break
            end
        end
        
        
    else
        
        while abs(f) > 1e-10
            
            % Transendental Equation
            [f,dfda,~]=transendentalfunc(a,N,alpha,beta,constants,mode);
            
            % Newton's method
            a = a - f/dfda;
            a = real(a);
            [alpha,beta,~,~,a_m,~,t_mN,tF]=Get_alpha_beta(a,s,c,theta,N,t_des,'lower',constants,mode);
            
            % Counter
            itr = itr + 1;
            
            % Max iterations
            if itr > maxiter
                disp('Hitting Namx max iter')
                break
                
            end
            
        end
        
    end
    
    % Clean up (miniscule img part)
    if tF < t_des
        Min_energy_Min_Time_Sols=vertcat(Min_energy_Min_Time_Sols,[N,a_m,t_mN,a,tF,t_des]);  %[N,a_m,t_mN,a_min,t_min,t_des]
    end
    
end

if Min_energy_Min_Time_Sols(end,5)>Min_energy_Min_Time_Sols(end,6)
    keyboard
end

output.Nmax=Nmax;
output.Min_energy_Min_Time_Sols=Min_energy_Min_Time_Sols;
output.parabolic=0;
end


function [f,dfda]=solve_for_amin(a,N,alpha,beta,constants,mode)
[f,dfda,~]=transendentalfunc(a,N,alpha,beta,constants,mode);

end