function Asol=Get_a(t_des,Nmax_Output,lambertparams,constants,mode)
% Find the a where the line t_des intersects the curves for all N
% The current implementation uses bisection method. Newtons method has slow
% convergence at a_m. Next version will have a polynomial/spline
% approximation to have a better initial condition.

% Asol is a matrix that saves the N, the corresponding a, and the branch.



theta=lambertparams.theta;
s=lambertparams.s;
c=lambertparams.c;

Nmax=Nmax_Output.Nmax;
Min_energy_Min_Time_Sols=Nmax_Output.Min_energy_Min_Time_Sols; %[N,a_m,t_mN,a_min,t_min]


% 3rd col is upper or lower branch
Asol=zeros(2*Nmax+1,4);% [N,a,branch,theta]
% upper : 1
% lower : -1

% Loop through all solution cases
ncnt=1;

%% solve for N==0

N=0;
MengyMtm=Min_energy_Min_Time_Sols(Min_energy_Min_Time_Sols(:,1)==N,:); %[N,a_m,t_mN,a_min,t_min]
a_m=MengyMtm(2);
t_mN=MengyMtm(3);
a_min=MengyMtm(4);
t_min=MengyMtm(5);

if t_des>=lambertparams.t_m
    a_tFsmall=a_m;
    a_tFlarge=a_m;
    % first march forward and get interval that has t_des
    while(1)
        a_tFsmall =a_tFlarge;
        a_tFlarge=a_tFlarge+0.2;
        [~,~,~,~,~,~,~,tF]=Get_alpha_beta(a_tFlarge,s,c,theta,N,t_des,'upper',constants,mode) ;
        if tF>t_des
            break
        end

    end
    a=bisection_method_for_a(t_des,N,a_tFsmall,a_tFlarge,'upper',Nmax_Output,lambertparams,constants,mode);
    Asol(ncnt,:)=[0,a,1,theta];
else
    a_tFsmall=a_m;
    a_tFlarge=a_m;
    % march downward till you enclose t_Des
    while(1)
        a_tFlarge=a_tFsmall;
        a_tFsmall=a_tFsmall+0.2;
        [~,~,~,~,~,~,~,tF]=Get_alpha_beta(a_tFsmall,s,c,theta,N,t_des,'lower',constants,mode) ;
        if tF<t_des
            break
        end

    end
    a=bisection_method_for_a(t_des,N,a_tFsmall,a_tFlarge,'lower',Nmax_Output,lambertparams,constants,mode);
    Asol(ncnt,:)=[0,a,-1,theta];
end



ncnt=ncnt+1;

%% solve for N>0
for N=1:Nmax
   
    % Determine starting "a" when t_min < t_des < t_m
    
    MengyMtm=Min_energy_Min_Time_Sols(Min_energy_Min_Time_Sols(:,1)==N,:); %[N,a_m,t_mN,a_min,t_min]
    a_m=MengyMtm(2);
    t_mN=MengyMtm(3);
    a_min=MengyMtm(4);
    t_min=MengyMtm(5);
    if t_des >t_mN 
        
        % upper branch ----------
        a_tFsmall=a_m;
        a_tFlarge=a_m;
        while(1)
            a_tFsmall=a_tFlarge;
            a_tFlarge=a_tFlarge+0.2;
            [~,~,~,~,~,~,~,tF]=Get_alpha_beta(a_tFlarge,s,c,theta,N,t_des,'upper',constants,mode) ;
            if tF>t_des
                break
            end
                
        end
        a=bisection_method_for_a(t_des,N,a_tFsmall,a_tFlarge,'upper',Nmax_Output,lambertparams,constants,mode);       
        Asol(ncnt,:)=[N,a,1,theta];
        
        % lower branch----------------
        a_tFsmall=a_min;
        a_tFlarge=a_min;
        while(1)
            a_tFsmall=a_tFlarge;
            a_tFlarge=a_tFlarge+0.2;
            [~,~,~,~,~,~,~,tF]=Get_alpha_beta(a_tFlarge,s,c,theta,N,t_des,'lower',constants,mode) ;
            if tF>t_des
                break
            end
                
        end
        a=bisection_method_for_a(t_des,N,a_tFsmall,a_tFlarge,'lower',Nmax_Output,lambertparams,constants,mode);       
        Asol(ncnt+1,:)=[N,a,-1,theta];
    
    else
        % left of a_min
        a_tFlarge=a_m;
        a_tFsmall=a_min;
        a=bisection_method_for_a(t_des,N,a_tFsmall,a_tFlarge,'lower',Nmax_Output,lambertparams,constants,mode);
        Asol(ncnt,:)=[N,a,-1,theta];

        % right of a_min
        a_tFsmall=a_min;
        a_tFlarge=a_min;
        while(1)
            a_tFsmall=a_tFlarge;
            a_tFlarge=a_tFlarge+0.2;
            [~,~,~,~,~,~,~,tF]=Get_alpha_beta(a_tFlarge,s,c,theta,N,t_des,'lower',constants,mode) ;
            if tF>t_des
                break
            end
                
        end
        a=bisection_method_for_a(t_des,N,a_tFsmall,a_tFlarge,'lower',Nmax_Output,lambertparams,constants,mode);
        Asol(ncnt+1,:)=[N,a,-1,theta];
        
    end

    ncnt=ncnt+2;
    
end

%         a0 = 0.95*a_m+0.05*a_min;
%         a0 = 1.05*a_min;

end