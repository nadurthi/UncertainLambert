
function perturbed_EFMfilteredsummary=Get_filtered_perturbed_EFM_summary(perturbed_efmoutput,filters,constants,config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% @Input:
% perturbed_EFMfilteredsummary: the efmoutput generated by
% EFMcompute_perturbed
% filters: The required filters. The filters should contain:
            % filters.a_min
            % filters.a_max
            % filters.delVD_min
            % filters.delVD_max
            % filters.Nmin_bnd
            % filters.Nmax_bnd
            % filters.min_alt
            % filters.max_alt

            % binary (0/1) : 0-> Disable, 1-> Enable filters
            % filters.check_escapevel  
            % filters.check_a_bnds
            % filters.check_delVD_bnds
            % filters.check_N_bnds
            % filters.check_alt_bnds

% @Output: 
%     perturbed_EFMfilteredsummary: The filtered EFM, summarized into
%                           - Number of solutions
%                           - min delV total
%                           - min delV departure
%                           - min delV arrival

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N_deptime=length(perturbed_efmoutput.inputconfig.DepTime);
N_timeflights=length(perturbed_efmoutput.inputconfig.TimeFlights);


[DT,TF]=meshgrid(perturbed_efmoutput.inputconfig.DepTime,perturbed_efmoutput.inputconfig.TimeFlights);

DDV_min=zeros(size(DT));
DDV_max=zeros(size(DT));

DDVDep_min=zeros(size(DT));
DDVDep_max=zeros(size(DT));

DDVArr_min=zeros(size(DT));
DDVArr_max=zeros(size(DT));

DNsol_all=zeros(size(DT));
DNsol_feas=zeros(size(DT));

DDV_min_pert=zeros(size(DT));
DDV_max_pert=zeros(size(DT));

DDVDep_min_pert=zeros(size(DT));
DDVDep_max_pert=zeros(size(DT));

DDVArr_min_pert=zeros(size(DT));
DDVArr_max_pert=zeros(size(DT));

DNsol_all_pert=zeros(size(DT));
DNsol_feas_pert=zeros(size(DT));

for depstarttime_ind = 1:N_deptime        % EFM x-axis (departure time)
    
    for timeflight_ind = 1:N_timeflights
        Asol=perturbed_efmoutput.SimSol_Asol{depstarttime_ind,timeflight_ind};
        Nsols_all=size(Asol,1);
        
        transferorbits=perturbed_efmoutput.SimSol_transferorbits_kepler{depstarttime_ind,timeflight_ind};
        transferorbits_pert=perturbed_efmoutput.SimSol_transferorbits{depstarttime_ind,timeflight_ind};
        
        delV=zeros(Nsols_all,1);
        delVD=zeros(Nsols_all,1);
        delVA=zeros(Nsols_all,1);
        
        delV_pert=zeros(Nsols_all,1);
        delVD_pert=zeros(Nsols_all,1);
        delVA_pert=zeros(Nsols_all,1);
        
        for i=1:Nsols_all
%             keyboard
            transorbit=convert_transferorbit_mat2struct(transferorbits(i,:));
            transorbit_pert=convert_transferorbit_mat2struct_perturbed(transferorbits_pert(i,:));
            
            [feas,~]=SolnFilter(transorbit,filters,constants,config);
            [feas_pert,~]=SolnFilter_pert(transorbit_pert,filters,constants,config);
            
            if i==1
                Feas=ones(Nsols_all,length(feas));
                Feas_pert=ones(Nsols_all,length(feas_pert));
            end
            Feas(i,:)=feas;
            Feas_pert(i,:)=feas_pert;
            
            delV(i)=transorbit.delV;
            delVD(i)=transorbit.delVD;
            delVA(i)=transorbit.delVA;
            
            delV_pert(i)=transorbit_pert.delV;
            delVD_pert(i)=transorbit_pert.delVD;
            delVA_pert(i)=transorbit_pert.delVA;
            
            
        end

        if isempty(Asol)==1
            DNsol_all(timeflight_ind,depstarttime_ind)=NaN;
            DNsol_feas(timeflight_ind,depstarttime_ind)=NaN;
            
            DDV_min(timeflight_ind,depstarttime_ind)=NaN;
            DDV_max(timeflight_ind,depstarttime_ind)=NaN;
            
            DDVDep_min(timeflight_ind,depstarttime_ind)=NaN;
            DDVDep_max(timeflight_ind,depstarttime_ind)=NaN;

            DDVArr_min(timeflight_ind,depstarttime_ind)=NaN;
            DDVArr_max(timeflight_ind,depstarttime_ind)=NaN;
            
            %
            DNsol_all_pert(timeflight_ind,depstarttime_ind)=NaN;
            DNsol_feas_pert(timeflight_ind,depstarttime_ind)=NaN;
            
            DDV_min_pert(timeflight_ind,depstarttime_ind)=NaN;
            DDV_max_pert(timeflight_ind,depstarttime_ind)=NaN;
            
            DDVDep_min_pert(timeflight_ind,depstarttime_ind)=NaN;
            DDVDep_max_pert(timeflight_ind,depstarttime_ind)=NaN;

            DDVArr_min_pert(timeflight_ind,depstarttime_ind)=NaN;
            DDVArr_max_pert(timeflight_ind,depstarttime_ind)=NaN;
            
            continue
        end
        

        feas=prod(Feas,2);
        feas_pert=prod(Feas_pert,2);
        
        DNsol_all(timeflight_ind,depstarttime_ind)=Nsols_all;
        DNsol_all_pert(timeflight_ind,depstarttime_ind)=Nsols_all;
        
        % take only feasible orbits
        Nsols_filtered=sum(feas);
        delV=delV(feas==1);
        delVD=delVD(feas==1);
        delVA=delVA(feas==1);
        
        
        if Nsols_filtered==0 % if empty i.e. no feasible solutions
            DNsol_feas(timeflight_ind,depstarttime_ind)=NaN;
            DDV_min(timeflight_ind,depstarttime_ind)=NaN;
            DDV_max(timeflight_ind,depstarttime_ind)=NaN;
            
            DDVDep_min(timeflight_ind,depstarttime_ind)=NaN;
            DDVDep_max(timeflight_ind,depstarttime_ind)=NaN;

            DDVArr_min(timeflight_ind,depstarttime_ind)=NaN;
            DDVArr_max(timeflight_ind,depstarttime_ind)=NaN;
        else
            DNsol_feas(timeflight_ind,depstarttime_ind)=Nsols_filtered;
            DDV_min(timeflight_ind,depstarttime_ind)=min( delV )*constants.normV2trueV;
            DDV_max(timeflight_ind,depstarttime_ind)=max( delV )*constants.normV2trueV;
            
            DDVDep_min(timeflight_ind,depstarttime_ind)=min( delVD )*constants.normV2trueV;
            DDVDep_max(timeflight_ind,depstarttime_ind)=max( delVD )*constants.normV2trueV;

            DDVArr_min(timeflight_ind,depstarttime_ind)=min( delVA )*constants.normV2trueV;
            DDVArr_max(timeflight_ind,depstarttime_ind)=max( delVA )*constants.normV2trueV;
        end
        
        
        % take only feasible orbits for perturbations
        delV_pert=delV_pert(feas_pert==1);
        delVD_pert=delVD_pert(feas_pert==1);
        delVA_pert=delVA_pert(feas_pert==1);
        
        delV_pert=delV_pert(~isnan(delV_pert));
        delVD_pert=delVD_pert(~isnan(delV_pert));
        delVA_pert=delVA_pert(~isnan(delV_pert));
        Nsols_filtered=length(delV_pert);
        
        if Nsols_filtered==0 % if empty i.e. no feasible solutions
            DNsol_feas_pert(timeflight_ind,depstarttime_ind)=NaN;
            DDV_min_pert(timeflight_ind,depstarttime_ind)=NaN;
            DDV_max_pert(timeflight_ind,depstarttime_ind)=NaN;
            
            DDVDep_min_pert(timeflight_ind,depstarttime_ind)=NaN;
            DDVDep_max_pert(timeflight_ind,depstarttime_ind)=NaN;

            DDVArr_min_pert(timeflight_ind,depstarttime_ind)=NaN;
            DDVArr_max_pert(timeflight_ind,depstarttime_ind)=NaN;
        else
            DNsol_feas_pert(timeflight_ind,depstarttime_ind)=Nsols_filtered;
            DDV_min_pert(timeflight_ind,depstarttime_ind)=min( delV_pert )*constants.normV2trueV;
            DDV_max_pert(timeflight_ind,depstarttime_ind)=max( delV_pert )*constants.normV2trueV;
            
            DDVDep_min_pert(timeflight_ind,depstarttime_ind)=min( delVD_pert )*constants.normV2trueV;
            DDVDep_max_pert(timeflight_ind,depstarttime_ind)=max( delVD_pert )*constants.normV2trueV;

            DDVArr_min_pert(timeflight_ind,depstarttime_ind)=min( delVA_pert )*constants.normV2trueV;
            DDVArr_max_pert(timeflight_ind,depstarttime_ind)=max( delVA_pert )*constants.normV2trueV;
        end
        
        
    end
end

perturbed_EFMfilteredsummary.DepTime=perturbed_efmoutput.inputconfig.DepTime;
perturbed_EFMfilteredsummary.TimeFlights=perturbed_efmoutput.inputconfig.TimeFlights;

perturbed_EFMfilteredsummary.DT=DT;
perturbed_EFMfilteredsummary.TF=TF;

perturbed_EFMfilteredsummary.DDV_min=DDV_min;
perturbed_EFMfilteredsummary.DDV_max=DDV_max;

perturbed_EFMfilteredsummary.DDVDep_min=DDVDep_min;
perturbed_EFMfilteredsummary.DDVDep_max=DDVDep_max;

perturbed_EFMfilteredsummary.DDVArr_min=DDVArr_min;
perturbed_EFMfilteredsummary.DDVArr_max=DDVArr_max;

perturbed_EFMfilteredsummary.DNsol_all=DNsol_all;
perturbed_EFMfilteredsummary.DNsol_feas=DNsol_feas;

%
perturbed_EFMfilteredsummary.DDV_min_pert=DDV_min_pert;
perturbed_EFMfilteredsummary.DDV_max_pert=DDV_max_pert;

perturbed_EFMfilteredsummary.DDVDep_min_pert=DDVDep_min_pert;
perturbed_EFMfilteredsummary.DDVDep_max_pert=DDVDep_max_pert;

perturbed_EFMfilteredsummary.DDVArr_min_pert=DDVArr_min_pert;
perturbed_EFMfilteredsummary.DDVArr_max_pert=DDVArr_max_pert;

perturbed_EFMfilteredsummary.DNsol_all_pert=DNsol_all_pert;
perturbed_EFMfilteredsummary.DNsol_feas_pert=DNsol_feas_pert;



