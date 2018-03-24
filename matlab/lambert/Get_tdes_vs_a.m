function X=Get_tdes_vs_a(N_max,lambertparams,constants)
% @brief : This function generates the cureves of time of flight vs
% semi-major axis a for given set of lambert parameters.

% @Input : N_max : the maximum number of revolutions (also a the same as number of curve+1)
%          lambertparams: contains r1,r2,s,c,theta and a_m( the minimum possible semi-major axis a)
% @Output: X: cell containing each curve


s=lambertparams.s;
c=lambertparams.c;
theta=lambertparams.theta;



Avec=[linspace(lambertparams.a_m,lambertparams.a_m*1.01,1000),linspace(lambertparams.a_m*1.01,lambertparams.a_m*1.2,1000),linspace(lambertparams.a_m*1.2,20,1000)];
X=cell(1,1);

for N=0:N_max
    P=[];
    % upper branch
    for a = Avec(end:-1:1)
        [alpha,beta,alpha_m,beta_m,a_m,t_m,t_mN,tF]=Get_alpha_beta(a,s,c,theta,N,0,'upper',constants);
        P=vertcat(P,[a,tF]);
    end
    X{N+1}=P;
    
    P=[];
    % upper branch
    for a = Avec
        [alpha,beta,alpha_m,beta_m,a_m,t_m,t_mN,tF]=Get_alpha_beta(a,s,c,theta,N,0,'lower',constants);
        P=vertcat(P,[a,tF]);
    end
    X{N+1}=vertcat(X{N+1},P);
    
    
end
