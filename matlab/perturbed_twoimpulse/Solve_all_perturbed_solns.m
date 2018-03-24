function pert_trajs=Solve_all_perturbed_solns(t1_true,t2_true,kepler_trajs,Satmodel,perturb_models,integrator,solver,solverprops,constants,config)
% This function call the user speicified solver, using the two-body
% transfer as a hot start.
% t1_true,t2_true are the start and end Ephemeris time : non-canonical
% units of time

% kepler_trajs: all the two-body two impulse solutions
% integrator: integrator to be used
% solver: solver to be used
% solverprops: solverprops has all the parameters required by the solver

for i=1:1:size(kepler_trajs,1)
    [i,size(kepler_trajs,1)]
    transorbit=convert_transferorbit_mat2struct(kepler_trajs(i,:));
    
    r1=transorbit.r1;
    r2=transorbit.r2;
    v1_0=transorbit.v1;
    
    
    %     keyboard
    
    
    normacc=true;
    propagator=@(Tvec,x0,epsilon)integrator(Tvec,x0,epsilon,perturb_models,normacc,constants,Satmodel,config);
    solverprops.Nt=100;
    NN=transorbit.N;
    
    h=cross(r1,v1_0);
    evec=cross(v1_0,h)/1-r1/norm(r1);
    e=norm(evec);
    disp(['NN,e = ',num2str(NN),' , ',num2str(e)])
    try
        [r1,rf2,v1,vf2,traj_t,traj_x]=solver(r1,r2,v1_0,t1_true*constants.trueT2normT,t2_true*constants.trueT2normT,propagator,solverprops);
        %         keyboard
        
    catch
        rnd=floor(10000*abs(randn(1)));
        save(strcat('Failed Case : ',datestr(now(),'dd-mmm-yyyy HH:MM:SS:FFF'),',',num2str(rnd)))
        disp(['User specified solver failed: ',' #revs = ',num2str(NN),' , eccentricity = ',num2str(e)])
        r1=r1;
        rf2=r2;
        v1=v1_0*NaN;
        vf2=v1_0*NaN;
        traj_t=linspace(t1_true*constants.trueT2normT,t2_true*constants.trueT2normT,solverprops.Nt);
        traj_x=ones(solverprops.Nt,6)*NaN;
        
    end
    
    
    
    pert_transorbit.r1=r1;
    pert_transorbit.r2=rf2(:)';
    pert_transorbit.v1=v1;
    pert_transorbit.v2=vf2(:)';
    pert_transorbit.timevec=traj_t;
    pert_transorbit.traj=traj_x;
    
    if any(isnan(pert_transorbit.v1))==0
        if config.useSPICE==true
            pert_transorbit.true_Altitude=get_sat_earth_altitude(pert_transorbit.timevec*constants.normT2trueT,pert_transorbit.traj(:,1:3)*constants.normX2trueX,constants);
        else
            pert_transorbit.true_Altitude=sqrt(sum(( pert_transorbit.traj(:,1:3)*constants.normX2trueX ).^2,2))-constants.Re;
        end
    else
        pert_transorbit.true_Altitude=pert_transorbit.timevec*NaN;
    end
    
    pert_transorbit.N= transorbit.N;
    
    pert_transorbit.v0A= transorbit.v0A;
    pert_transorbit.v0B= transorbit.v0B;
    
    pert_transorbit.delV = norm(abs(pert_transorbit.v1-pert_transorbit.v0A)) + norm(abs(pert_transorbit.v2-pert_transorbit.v0B)); % Delta V cost #changed norm to 1 norm
    pert_transorbit.delVD = norm(abs(pert_transorbit.v1-pert_transorbit.v0A)); % Departure
    pert_transorbit.delVA = norm(abs(pert_transorbit.v2-pert_transorbit.v0B)); % Arrival
    
    M=convert_transferorbit_struct2mat_perturbed(pert_transorbit);
    if i==1
        pert_trajs=zeros(size(kepler_trajs,1),length(M));
    end
    pert_trajs(i,:)=M;
    

end