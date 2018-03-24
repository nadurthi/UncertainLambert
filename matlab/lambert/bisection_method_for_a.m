function a=bisection_method_for_a(t_des,N,a_tFsmall,a_tFlarge,branch,Nmax_Output,lambertparams,constants,config)
% @brief : Bisection method to solve for a, given a t_des
% @ Input: t_des - time of flight
%          N - the number of revolutions
%          a_tFsmall - the semi-major axis a for which tf<t_des
%          a_tFsmall - the semi-major axis a for which tf>t_des
%          branch - lower or upper branch
%          lambertparams - lambertparams is set of lambert parameters

% @Output : return the a semi-major axis for a given t_des,N and brnach

a_m=lambertparams.a_m;
t_m=lambertparams.t_m;
s=lambertparams.s;
c=lambertparams.c;
theta=lambertparams.theta;


a1=a_tFsmall;
a2=a_tFlarge;
itr=1;
while(1)
    amid=(a1+a2)/2;
    if strcmp(branch,'upper')
        [~,~,~,~,~,~,~,tF]=Get_alpha_beta(amid,s,c,theta,N,t_des,'upper',constants) ;
    elseif strcmp(branch,'lower')
        [~,~,~,~,~,~,~,tF]=Get_alpha_beta(amid,s,c,theta,N,t_des,'lower',constants) ;
    else
        error('Only use upper or lower')
    end
    
    if abs(tF-t_des)<=1e-9
        break
    end
    if tF<t_des
        a1=amid;
    else
        a2=amid;
    end
    itr=itr+1;

    if itr>150
        branch
        [abs(tF-t_des),N,amid,a_tFsmall,a_tFlarge]
        error('bisection itr>150')
    end
end
a=(a1+a2)/2;
a=real(a);
