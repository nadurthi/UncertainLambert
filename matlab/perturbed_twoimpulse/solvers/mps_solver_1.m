function [r1,rf2,v1_0,vf2,traj_t,traj_x]=mps_solver_1(r1,r2,v1_0,t1,t2,propagator,solverprops)
% use method of particular solutions to get initial velocity
% v1_0 is the initial velocity estimate
% propagator is called as [t,x]=propagator(T,x0,epsilon)

r1=r1(:);
r2=r2(:);
v1_0=v1_0(:);

delVvecs=solverprops.mps_1.get_delVs(solverprops.mps_1.Np,1); % V=[... delv1 ...;
%                                         ... delv2 ...;...]

err=1000;



while (err>solverprops.mps_1.tol)
    delVvecs=delVvecs*exp((1-iter)/(0.3*solverprops.mps_1.maxiter));
    [t,xref]=propagator([t1,t2],[r1;v1_0],1);
    
    %     Xf=cell(1,size(delVvecs,1));
    A=zeros(3,size(delVvecs,1));
    DV=zeros(3,size(delVvecs,1));
    for i=1:size(delVvecs,1)
        v=v1_0+delVvecs(i,:)';
        [t,x]=propagator([t1,t2],[r1;v],1);
        
        DV(:,i)=delVvecs(i,:)';
        A(:,i)=x(end,1:3)'-xref(end,1:3)';
    end
    alphas=pinv(A)*(r2-xref(end,1:3)');
    v1_0=v1_0+DV*alphas;
    
    [t,x]=propagator([t1,t2],[r1,v1_0],1);
    rf2=x(end,1:3)';
    vf2=x(end,4:6)';
    err=norm(r2-rf2);
    %             [epsilon,iter,err]
    iter=iter+1;
    
    if iter>solverprops.mps_1.maxiter
        error('MPS failed: reached maxiter')
    end
    
end

[traj_t,traj_x]=propagator(linspace(t1,t2,solverprops.Nt),[r1,v1_0],1);
rf2=traj_x(end,1:3)';
vf2=traj_x(end,4:6)';
