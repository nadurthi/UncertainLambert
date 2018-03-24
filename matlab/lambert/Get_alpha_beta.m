function [alpha,beta,alpha_m,beta_m,a_m,t_m,t_mN,tF]=Get_alpha_beta(a,s,c,theta,N,t_des,branch,constants,config)
% alpha_m, beta_m : minimum energy alpha and beta
% a_m,t_m : minimum energy transfer a and time of transfer
% t_mN: minimum energy time of transfer for given N
% tF: time of flight for given a, branch and N.

% branch={'upper' or 'lower' or NaN}
% if branch = NaN then t_des is used to decide which branch: this is not
% used anywhere ... might be deprecated later

a_m=s/2;
alpha_m0=pi;
beta_m0  = 2*asin(sqrt((s-c)/(s)));
if theta >= 0 && theta <= pi
    beta_m = beta_m0;                      % Eq. 5.31a
elseif theta >= pi && theta < 2*pi
    beta_m = -beta_m0;                     % Eq. 5.31b
end
alpha_m=alpha_m0;

t_m=1/sqrt(constants.muCan)*(s/2)^(3/2)*( pi-beta_m+sin(beta_m) );
t_mN=1/sqrt(constants.muCan)*(s/2)^(3/2)*( (2*N+1)*pi-beta_m+sin(beta_m) );

alpha_0 = 2*asin(sqrt(s/(2*a)));        % Eq. 5.26
beta_0  = 2*asin(sqrt((s-c)/(2*a)));    % Eq. 5.27

% Beta
if theta >= 0 && theta < pi
    beta = beta_0;                      % Eq. 5.31a
elseif theta >= pi && theta < 2*pi
    beta = -beta_0;                     % Eq. 5.31b
end

% Alpha

if any(isnan(branch)==1)
    if t_des <= t_mN
        alpha = alpha_0;                    % Eq. 5.32a
    elseif t_des > t_mN
        alpha = 2*pi - alpha_0;             % Eq. 5.32b
    end
    
elseif strcmp(branch,'lower')
    alpha = alpha_0;

elseif strcmp(branch,'upper')
    alpha = 2*pi - alpha_0;
    
else
    error('Option not valid')
end

alpha  = real(alpha);
beta   = real(beta);


tF = sqrt((a^3)/constants.muCan)*(2*N*pi + alpha - beta - (sin(alpha) - sin(beta))); % Eq. 5.41

 