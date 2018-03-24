function [r1,rf2,v1,vf2,traj_t,traj_x]=shooting_solver(r1,r2,v1_0,t1,t2,propagator,solverprops)
% use method of particular solutions to get initial velocity
% v1_0 is the initial velocity estimate
% propagator is called as [t,x]=propagator(T,x0,epsilon)

options = optimoptions(@fsolve,'FunctionTolerance',1e-6,'StepTolerance',1e-6,'Display','iter','MaxIterations',25,'MaxFunctionEvaluations',10000); %,'Algorithm','levenberg-marquardt'

[v1,fval,exitflag,output] =fsolve(@(v)shoot_v1(v,r1,r2,t1,t2,@(t,x)propagator(t,x,1)),v1_0,options);

if exitflag<=0
    error('Max iteration reached for shooting')
end


if exitflag<=0
    options = optimoptions(@fsolve,'FunctionTolerance',1e-6,'StepTolerance',1e-6,'Display','off','MaxIterations',25,'MaxFunctionEvaluations',10000); %,'Algorithm','levenberg-marquardt'
    disp('Switching to homotopy approach')
    % do homotopy approach
%     epsilon_vec=[linspace(0.1,0.65,5),linspace(0.7,1,10)];
    epsilon_vec=linspace(0.1,1,10);
    v1=v1_0;
    for epsilon=epsilon_vec
%         epsilon
        [v1,fval,exitflag,output] =fsolve(@(v)shoot_v1(v,r1,r2,t1,t2,@(t,x)propagator(t,x,epsilon)),v1,options);
        if exitflag<=0
            error('Error homotopy shooting fsolve failed')
        end
    end
    
    
end


[traj_t,traj_x]=propagator(linspace(t1,t2,solverprops.Nt),[r1,v1],1);
rf2=traj_x(end,1:3)';
vf2=traj_x(end,4:6)';
end

function err=shoot_v1(v,r1,r2,t1,t2,propagator)
[t,x]=propagator([t1,t2],[r1,v]);
err=(r2-x(end,1:3));

end




