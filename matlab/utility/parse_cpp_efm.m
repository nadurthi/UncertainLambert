function efmoutput=parse_cpp_efm(efmsize,efmoutput,L)

% efmoutput.ArrOrbit_normalized_cell=ArrOrbit_normalized_cell;
% efmoutput.inputconfig=inputconfig;
% efmoutput.SimSol_Nmax=SimSol_Nmax;
% efmoutput.SimSol_DeltaVelocities=SimSol_DeltaVelocities;
% efmoutput.SimSol_Asol=SimSol_Asol;
%
% efmoutput.inputconfig=inputconfig;
% efmoutput.SimSol_transferorbits=SimSol_transferorbits;


depstarttime_ind=1;
timeflight_ind=1;
for ind=1:efmsize
    [ind,depstarttime_ind,timeflight_ind]
    %     {depstarttime_ind,timeflight_ind}
    %     SimSol_Nmax
    %     SimSol_DeltaVelocities
    %     SimSol_Asol
    %     SimSol_transferorbits
    
    t_des = str2double(L{ind}.t_des);
    efmoutput.inputconfig.DepOrbit_normalized{depstarttime_ind} = convertstrrow2doublevector(L{ind}.xA);
    efmoutput.inputconfig.ArrOrbit_normalized_cell{depstarttime_ind,timeflight_ind} = {convertstrrow2doublevector(L{ind}.xB),t_des};
    
    %%
    %      LP.mag_r1 << ',' << LP.mag_r2 << ','<< LP.c << ','<< LP.s << ','<< LP.a_m << ','<< LP.theta << ','<< LP.t_p << ','<< LP.t_m << "\n";
    LP1 = convertstrrow2doublevector(L{ind}.LP{1});
    theta1 = LP1(6);
    NN= length(L{ind}.LCPV{1});
    Nmax = (NN-1)/2;
    Snmax=[];
    Snmax.Nmax = Nmax;
    if Snmax.Nmax<0
        Snmax.parabolic=1;
    else
        Snmax.parabolic=0;
    end
    Snmax.Min_energy_Min_Time_Sols=zeros(NN,6); %[N,a_m,t_mN,a_min,t_min,t_des]
    for p = 1:length(L{ind}.LCPV{1})
        
        M=convertstrrow2doublevector(L{ind}.LCPV{1}{p}) ;%LCP.a_m 		<< ',' << LCP.t_mN 	<< ',' << LCP.a_min << ','<< LCP.t_min 	<< ',' << LCP.N
        try
            Snmax.Min_energy_Min_Time_Sols(p,:) = [M(5),M(1),M(2),M(3),M(4),t_des];
        catch
            keyboard
        end
    end
    efmoutput.inputconfig.SimSol_Nmax{depstarttime_ind,timeflight_ind}={Snmax,theta1 };
    
    LP2 = convertstrrow2doublevector(L{ind}.LP{2});
    theta2 = LP2(6);
    NN= length(L{ind}.LCPV{2});
    Nmax = (NN-1)/2;
    Snmax=[];
    Snmax.Nmax = Nmax;
    if Snmax.Nmax<0
        Snmax.parabolic=1;
    else
        Snmax.parabolic=0;
    end
    for p = 1:length(L{ind}.LCPV{2})
        M=convertstrrow2doublevector(L{ind}.LCPV{2}{p}) ;%LCP.a_m 		<< ',' << LCP.t_mN 	<< ',' << LCP.a_min << ','<< LCP.t_min 	<< ',' << LCP.N
        Snmax.Min_energy_Min_Time_Sols(p,:) = [M(5),M(1),M(2),M(3),M(4),t_des];
    end
    
    efmoutput.inputconfig.SimSol_Nmax{depstarttime_ind,timeflight_ind}{2,1}=Snmax;
    efmoutput.inputconfig.SimSol_Nmax{depstarttime_ind,timeflight_ind}{2,2}=theta2;
    
    %% Asol
    %     ss << LS.a 		<< ',' << LS.t_des 	<< ',' << LS.theta << ','<< LS.N 	<< "," << 1 << '\n';
    %     {depstarttime_ind,timeflight_ind}=Asol;%[N,a,branch,theta]
    
    if  efmoutput.inputconfig.SimSol_Nmax{depstarttime_ind,timeflight_ind}{1,1}.parabolic ==0
        NN= length(L{ind}.LSV{1});
        
        try
            Asol=zeros(NN,4); %[N,a,branch,theta]
            for p = 1:NN
                M=convertstrrow2doublevector(L{ind}.LSV{1}{p}) ;
                Asol(p,:) = [M(4),M(1),M(5), LP1(6)];
            end
            efmoutput.inputconfig.SimSol_Asol{depstarttime_ind,timeflight_ind}=Asol;
        catch
            keyboard
        end
    end
    if  efmoutput.inputconfig.SimSol_Nmax{depstarttime_ind,timeflight_ind}{2,1}.parabolic ==0
        NN= length(L{ind}.LSV{2});
        
        try
            Asol=zeros(NN,4); %[N,a,branch,theta]
            for p = 1:NN
                M=convertstrrow2doublevector(L{ind}.LSV{2}{p}) ;
                Asol(p,:) = [M(4),M(1),M(5), LP2(6)];
            end
            
            efmoutput.inputconfig.SimSol_Asol{depstarttime_ind,timeflight_ind}=vertcat(efmoutput.inputconfig.SimSol_Asol{depstarttime_ind,timeflight_ind},Asol);
        catch
            keyboard
        end
    end
    
    %% transfer orbit
    if  efmoutput.inputconfig.SimSol_Nmax{depstarttime_ind,timeflight_ind}{1,1}.parabolic ==0
        if  efmoutput.inputconfig.SimSol_Nmax{depstarttime_ind,timeflight_ind}{1,1}.parabolic ==0
            NN= length(L{ind}.LSV{1});
            
            
            for p = 1:NN
                transorbit=[];
                M=convertstrrow2doublevector(L{ind}.LTV{1}{p+1-1}) ;
                transorbit.r1=M(1:3);
                transorbit.v1=M(4:6);
                
                M=convertstrrow2doublevector(L{ind}.LTV{1}{p+2-1}) ;
                transorbit.r2=M(1:3);
                transorbit.v2=M(4:6);
                
                M=convertstrrow2doublevector(L{ind}.LTV{1}{p+3-1}) ;
                transorbit.v0A=M;
                
                M=convertstrrow2doublevector(L{ind}.LTV{1}{p+4-1}) ;
                transorbit.v0B=M;
                
                M=convertstrrow2doublevector(L{ind}.LTV{1}{p+5-1}) ;
                transorbit.delVD=M(1);
                transorbit.delVA=M(2);
                transorbit.delV=M(3);
                
                M=convertstrrow2doublevector(L{ind}.LTV{1}{p+6-1}) ;
                time=M;
                transorbit.normalized_transfer_timesteps=efmoutput.inputconfig.DepTime(depstarttime_ind)+time;
                
                M1=convertstrrow2doublevector(L{ind}.LTV{1}{p+7-1}) ;
                M2=convertstrrow2doublevector(L{ind}.LTV{1}{p+8-1}) ;
                transorbit.traj=[M1;M2];
                
                transorbit.t_des = t_des;
                transorbit.theta = theta1;
                try
                    efmoutput.inputconfig.SimSol_transferorbits{depstarttime_ind,timeflight_ind}{p}=transorbit;
                catch
                    keyboard
                end
            end
            
        end
    end
    if  efmoutput.inputconfig.SimSol_Nmax{depstarttime_ind,timeflight_ind}{2,1}.parabolic ==0
        Np=p;
        if  efmoutput.inputconfig.SimSol_Nmax{depstarttime_ind,timeflight_ind}{2,1}.parabolic ==0
            NN= length(L{ind}.LSV{2});
            
            
            for p = 1:NN
                transorbit=[];
                M=convertstrrow2doublevector(L{ind}.LTV{2}{p+1-1}) ;
                transorbit.r1=M(1:3);
                transorbit.v1=M(4:6);
                
                M=convertstrrow2doublevector(L{ind}.LTV{2}{p+2-1}) ;
                transorbit.r2=M(1:3);
                transorbit.v2=M(4:6);
                
                M=convertstrrow2doublevector(L{ind}.LTV{2}{p+3-1}) ;
                transorbit.v0A=M;
                
                M=convertstrrow2doublevector(L{ind}.LTV{2}{p+4-1}) ;
                transorbit.v0B=M;
                
                M=convertstrrow2doublevector(L{ind}.LTV{2}{p+5-1}) ;
                transorbit.delVD=M(1);
                transorbit.delVA=M(2);
                transorbit.delV=M(3);
                
                M=convertstrrow2doublevector(L{ind}.LTV{2}{p+6-1}) ;
                time=M;
                transorbit.normalized_transfer_timesteps=efmoutput.inputconfig.DepTime(depstarttime_ind)+time;
                
                M1=convertstrrow2doublevector(L{ind}.LTV{2}{p+7-1}) ;
                M2=convertstrrow2doublevector(L{ind}.LTV{2}{p+8-1}) ;
                transorbit.traj=[M1;M2];
                
                transorbit.t_des = t_des;
                transorbit.theta = theta1;
                
                efmoutput.inputconfig.SimSol_transferorbits{depstarttime_ind,timeflight_ind}{Np+p}=transorbit;
            end
            
        end
    end
    
    %%
    if timeflight_ind==efmoutput.inputconfig.TFgridlen
        depstarttime_ind=depstarttime_ind+1;
        timeflight_ind=0;
    end
    
    timeflight_ind=timeflight_ind+1;
    
end


