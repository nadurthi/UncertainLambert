function [t,x]=matlab_ode45_integrator(T,x0,epsilon,perturb_models,normacc,constants,Satmodel,config)
% T is vector with T(1)=t0 for given pair(t0,x0)
% x0 is initial state for given pair (t0,x0)
% epsilon is homotopy parameter for perturbed accelerations: 
%                                            epsilon =0 --> no perturbations
%                                            epsilon =1 --> full perturbations
% perturb_models is a cell of perturbation models.
% normacc=true --> normalize the accelerations (using canonical units)
% normacc=false --> integrate the ode in standard units km, km/s, km/s^2


x0=x0(:);

if epsilon>1 || epsilon<0
    error('epsilon should be between 0 and 1')
end

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
if normacc==true
    [t,x]=ode45(@(t,x)getxdot_normalized(t,x,epsilon,perturb_models,constants,Satmodel,config),T,x0,options);
else
    [t,x]=ode45(@(t,x)getxdot(t,x,epsilon,perturb_models,constants,Satmodel,config),T,x0,options);
end

end


function xdot=getxdot(t,x,epsilon,perturb_models,constants,Satmodel,config)

xyz_acc_pert = epsilon*get_total_perturbations(t,x,perturb_models,constants,Satmodel,config);

r=norm(x(1:3));

xdot(1,1)=x(4);
xdot(2,1)=x(5);
xdot(3,1)=x(6);
xdot(4,1)=-constants.mu*x(1)/r^3+xyz_acc_pert(1);
xdot(5,1)=-constants.mu*x(2)/r^3+xyz_acc_pert(2);
xdot(6,1)=-constants.mu*x(3)/r^3+xyz_acc_pert(3);

end

function xdot=getxdot_normalized(t,x,epsilon,perturb_models,constants,Satmodel,config)


y(1:3)=x(1:3)*constants.normX2trueX;
y(4:6)=x(4:6)*constants.normV2trueV;
xyz_acc_pert = epsilon*get_total_perturbations(t*constants.normT2trueT,y,perturb_models,constants,Satmodel,config);
xyz_acc_pert_normalized=xyz_acc_pert*constants.trueA2normA;


r=norm(x(1:3));


xdot(1,1)=x(4);
xdot(2,1)=x(5);
xdot(3,1)=x(6);
xdot(4,1)=-1*x(1)/r^3+xyz_acc_pert_normalized(1);
xdot(5,1)=-1*x(2)/r^3+xyz_acc_pert_normalized(2);
xdot(6,1)=-1*x(3)/r^3+xyz_acc_pert_normalized(3);

end

