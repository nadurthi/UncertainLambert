function perturbed_efmoutput=EFMcompute_perturbed(savefilename,Dep_orb_elem,dep_epoch,Arr_orb_elem,arr_epoch,depgridlen,TFgridlen,dept0,deptf,TFt0,TFtf,Satmodel,perturb_models,integrator,solver,solverprops,constants,config)


% Satmodel: Satellite model parameters
% perturb_models: cell of functions that return perturbed accelerations
% integrator: integrator function
% solver: solver function
% solverprops: the parameters for the solver

inputconfig.Dep_orb_elem=Dep_orb_elem;
inputconfig.Arr_orb_elem=Arr_orb_elem;

% Departure Orbit Elements
inputconfig.DepA    = Dep_orb_elem(1);             % km
inputconfig.DepE    = Dep_orb_elem(2);              % eccentricity
inputconfig.DepI    = Dep_orb_elem(3);          % Inclination (rad)
inputconfig.DepOm   = Dep_orb_elem(4);          % Right Ascension of Ascending Node (rad)
inputconfig.DepW    = Dep_orb_elem(5);          % Argument of Perigee (rad)
inputconfig.DepM0   = Dep_orb_elem(6);          % Mean Anomaly (rad)
[inputconfig.rD,inputconfig.vD] = elm2rv(inputconfig.DepA,inputconfig.DepE,inputconfig.DepI,inputconfig.DepOm,inputconfig.DepW,inputconfig.DepM0,0,constants.mu);
inputconfig.DepP    = 2*pi*sqrt(inputconfig.DepA^3 / constants.mu);
inputconfig.dep_epoch=dep_epoch;

% Arrival Orbit Elements
inputconfig.ArrA    = Arr_orb_elem(1);             % km
inputconfig.ArrE    = Arr_orb_elem(2);               % eccentricity
inputconfig.ArrI    = Arr_orb_elem(3);         % Inclination (rad)
inputconfig.ArrOm   = Arr_orb_elem(4);          % Right Ascension of Ascending Node (rad)
inputconfig.ArrW    = Arr_orb_elem(5);          % Argument of Perigee (rad)
inputconfig.ArrM0   = Arr_orb_elem(6);          % Mean Anomaly (rad)
[inputconfig.rA,inputconfig.vA] = elm2rv(inputconfig.ArrA,inputconfig.ArrE,inputconfig.ArrI,inputconfig.ArrOm,inputconfig.ArrW,inputconfig.ArrM0,0,constants.mu);
inputconfig.ArrP    = 2*pi*sqrt(inputconfig.ArrA^3 / constants.mu);
inputconfig.arr_epoch=arr_epoch;



%% Compute Time Vectors
inputconfig.dept0       = dept0;                            % Initial Departure Time
inputconfig.deptf      = deptf;                    % Final Departure Time


inputconfig.TFt0   = TFt0;     % min Time of Flight (y-axis)
inputconfig.TFtf   = TFtf;     %max Time of Flight

inputconfig.DepTime  = linspace(dept0,deptf,depgridlen);                  % Vector of Departure Times
inputconfig.TimeFlights = linspace(TFt0,TFtf,TFgridlen);                % Vector of Times of Flight

% Generate Departure & Arrival Orbits

propagator=@(T,x0)integrator(T,x0,1,perturb_models,false,constants,Satmodel,config);
[TimeFlights_normalized,DepOrbit_normalized,ArrOrbit_normalized_cell]=GetDepArr_normalized_orbit_positions_perturbed(inputconfig,propagator,constants,config);

% keyboard
%% EFM Generation
% save the solution for all departure times and time of flights

N_deptime=length(inputconfig.DepTime);
N_timeflights=length(inputconfig.TimeFlights);

SimSol_Nmax=cell(N_deptime,N_timeflights); %max # of solutions
SimSol_Asol=cell(N_deptime,N_timeflights); % the semi major axis a for each N

SimSol_transferorbits_kepler=cell(N_deptime,N_timeflights);
SimSol_transferorbits=cell(N_deptime,N_timeflights);
SimSol_DeltaVelocities=cell(N_deptime,N_timeflights); % total velocity change required: departure and arrival

tic


inputconfig.TimeFlights_normalized=TimeFlights_normalized;
inputconfig.DepOrbit_normalized=DepOrbit_normalized;
inputconfig.ArrOrbit_normalized_cell=ArrOrbit_normalized_cell;

%% convert everything to cells for parfor
DepOrbit_normalized_cell=cell(1,size(DepOrbit_normalized,1));
for depstarttime_ind = 1:N_deptime 
    DepOrbit_normalized_cell{depstarttime_ind}= DepOrbit_normalized(depstarttime_ind,:);
end

TU=constants.TU;
DepTime=inputconfig.DepTime;
spicepath=config.spicepath;
useSPICE=config.useSPICE;

%% Run in parallel
% now run for grid of efm
disp('Running for each combo')
parfor depstarttime_ind = 1:N_deptime        % EFM x-axis (departure time)
    if useSPICE
        load_spice_kernels(spicepath)
    end
    % Initial Position (Orbit A)
    r0A = DepOrbit_normalized_cell{depstarttime_ind}(1:3);
    v0A = DepOrbit_normalized_cell{depstarttime_ind}(4:6);
    deptime_normalized=DepTime(depstarttime_ind)/TU;
    
    % parfor computation
    for timeflight_ind = 1:N_timeflights    % EFM y-axis (time-of-flight)
%         [depstarttime_ind,timeflight_ind]
        disp(['On index : ',num2str(depstarttime_ind),' , ',num2str(timeflight_ind)])
        t_des = TimeFlights_normalized(timeflight_ind);
        
        
        X=ArrOrbit_normalized_cell{depstarttime_ind}{timeflight_ind};
        r0B = X(1:3);
        v0B = X(4:6);
        
        r1=norm(r0A);
        r2=norm(r0B);
        theta      = acos(r0A*r0B'/(r1*r2));

        % for a given r1,r2,theta ... get lambert problem parameters  
        lambertparams=GetLambert_parameters(r1,r2,theta,constants,config);

        % first find Nmax for the given t_des
        Nmax_Output=Get_Nmax(t_des,lambertparams,constants,config);
        SimSol_Nmax{depstarttime_ind,timeflight_ind}={theta,Nmax_Output};
        

        % then find 2*Nmax+1 a solutions for the given Nmax
        if Nmax_Output.parabolic ==0 && isfinite(Nmax_Output.Nmax) % Parabolic check
            
            %get the 2*Nmax+1 number of a, corresponding branch
            Asol=Get_a(t_des,Nmax_Output,lambertparams,constants,config);
            SimSol_Asol{depstarttime_ind,timeflight_ind}=Asol;
            
            % Solve for the trajectories
            kepler_trajs=Get_all_transferorbit_properties(Asol,r0A,r0B,v0A,v0B,deptime_normalized,t_des,lambertparams,constants,config);
            SimSol_transferorbits_kepler{depstarttime_ind,timeflight_ind}=kepler_trajs;
            
            % Solve for perturbed solutions
            t1_true=deptime_normalized*TU;
            t2_true=deptime_normalized*TU+t_des*TU;
            SimSol_transferorbits{depstarttime_ind,timeflight_ind}=Solve_all_perturbed_solns(t1_true,t2_true,kepler_trajs,Satmodel,perturb_models,integrator,solver,solverprops,constants,config);
        end
        
        
        
        theta =2*pi-theta; %is also a possibility, next version considers
        
        
        % for a given r1,r2,theta ... get all the parameters for the
        % lambert problem
        lambertparams=GetLambert_parameters(r1,r2,theta,constants,config);
        
        % first find Nmax for the given t_des
        Nmax_Output=Get_Nmax(t_des,lambertparams,constants,config);
        SimSol_Nmax{depstarttime_ind,timeflight_ind}{2,1}=theta;
        SimSol_Nmax{depstarttime_ind,timeflight_ind}{2,2}=Nmax_Output;
        
        % then find 2*Nmax+1 a solutions for the given Nmax
        if Nmax_Output.parabolic ==0 && isfinite(Nmax_Output.Nmax) % Parabolic check
            
            %get the 2*Nmax+1 number of a, corresponding branch
            Asol=Get_a(t_des,Nmax_Output,lambertparams,constants,config);
            SimSol_Asol{depstarttime_ind,timeflight_ind}=vertcat(SimSol_Asol{depstarttime_ind,timeflight_ind},Asol);
            
            % Solve for the trajectories
            kepler_trajs=Get_all_transferorbit_properties(Asol,r0A,r0B,v0A,v0B,deptime_normalized,t_des,lambertparams,constants,config);
            SimSol_transferorbits_kepler{depstarttime_ind,timeflight_ind}=[SimSol_transferorbits_kepler{depstarttime_ind,timeflight_ind};kepler_trajs];
            
            % Solve for perturbed solutions
            t1_true=deptime_normalized*TU;
            t2_true=deptime_normalized*TU+t_des*TU;
            pert_trajs=Solve_all_perturbed_solns(t1_true,t2_true,kepler_trajs,Satmodel,perturb_models,integrator,solver,solverprops,constants,config);
            SimSol_transferorbits{depstarttime_ind,timeflight_ind}=[SimSol_transferorbits{depstarttime_ind,timeflight_ind};pert_trajs];
        end
        
    end
    
    
end
toc

%% Saving

perturbed_efmoutput.ArrOrbit_normalized_cell=ArrOrbit_normalized_cell;
perturbed_efmoutput.inputconfig=inputconfig;
perturbed_efmoutput.SimSol_Nmax=SimSol_Nmax;
perturbed_efmoutput.SimSol_DeltaVelocities=SimSol_DeltaVelocities;
perturbed_efmoutput.SimSol_Asol=SimSol_Asol;

perturbed_efmoutput.inputconfig=inputconfig;
perturbed_efmoutput.SimSol_transferorbits_kepler=SimSol_transferorbits_kepler;
perturbed_efmoutput.SimSol_transferorbits=SimSol_transferorbits;

if isnan(savefilename)==0
    save(savefilename,'perturbed_efmoutput','-v7.3')
end