function output=GetLambert_parameters(r1,r2,theta,constants,mode)
% Get all the lambert constants for a fixed r1,r2,theta
% r1 and r2 are scalars

c       = sqrt(r1^2+r2^2-2*r1*r2*cos(theta));                          % Chord
s       = (r1 + r2 + c)/2;                        % Semiperimeter
a_m     = s/2;                                    % Minimum energy semimajor axis



        
t_p = sqrt(2)/3*(s^1.5 -sign(sin(theta))*(s-c)^1.5)/sqrt(constants.muCan);   % Eq. 5.34



beta_m0  = 2*asin(sqrt((s-c)/(s)));
if theta >= 0 && theta <= pi
    beta_m = beta_m0;                      % Eq. 5.31a
elseif theta >= pi && theta <= 2*pi
    beta_m = -beta_m0;                     % Eq. 5.31b
end


t_m=1/sqrt(constants.muCan)*(s/2)^(3/2)*( pi-beta_m+sin(beta_m) );

output.t_m=t_m;
output.r1=r1;
output.r2=r2;
output.c=c;
output.s=s;
output.a_m=a_m;
output.theta=theta;
output.t_p=t_p;

end