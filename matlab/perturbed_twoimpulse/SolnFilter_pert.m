function [feas,FeasCodes_index]=SolnFilter_pert(transorbit,filters,constants,mode)
% transorbit contains the transfer orbit trajectory and its properties in "normalized" form

% filters should contain these variables
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



% neg_a : col1
% max_a : col2
% earthcol : col3
% escapevel : col4
% delVfinite : col5
% transfer-arrival-match : col6
% a>=a_m : col7
% the a actually results in t_des : col8
% initial delv has to be less than 11km/s

FeasCodes_index.neg_a=1;
FeasCodes_index.a_isin_rng=2;
FeasCodes_index.earthcol=3;
FeasCodes_index.escapevel=4;
FeasCodes_index.delVfinite=5;
FeasCodes_index.transfer_arrival_match=6;
FeasCodes_index.a_gte_a_m=7;
FeasCodes_index.a_match_t_des=8;
FeasCodes_index.delVdep_max=9;
FeasCodes_index.N_isin_rng=10;
FeasCodes_index.delVarr_max=11;

% initialize the feasibility vector.
feas=zeros(1,11);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  check negative a
feas(1)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check_a_bnds
feas(2)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check Earth Collision
if filters.check_alt_bnds==1
%     Altitude=get_sat_earth_altitude(transorbit.normalized_transfer_timesteps*constants.TU,transorbit.traj(:,1:3)*constants.normX2trueX,constants);
%     [min(Altitude),max(Altitude)]
    if any(transorbit.true_Altitude < filters.min_alt) || any(transorbit.true_Altitude > filters.max_alt)
        feas(3)=0;
    else
        feas(3)=1;
    end
else
    feas(3)=1; % if not checking then take it to be feasible
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check EscVel
if filters.check_escapevel==1
    normV=sqrt(sum(transorbit.traj(:,4:6).^2,2))*constants.normV2trueV;
    if any(normV)>=11.2
        feas(4)=0;
    else
        feas(4)=1;
    end
else
    feas(4)=1;
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if delV is not finite
if isfinite(transorbit.delV)==0
    feas(5)=0;
else
    feas(5)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the transfer orbit acutally intersect the arrival orbit
if norm(transorbit.traj(1,1:3)-transorbit.r1(:)')+norm(transorbit.traj(end,1:3)-transorbit.r2(:)')>1e-5
    disp('failed transfer r1 and r2 do not match to TF orbit positions')
    feas(6)=0;
else
    feas(6)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cannot be less than a_m
feas(7)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if the timeofflight from a is acutally close to t_des
feas(8)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if departure delv is less than 11km/s
% remember delVD is normalized, so bring it to original km/s
if filters.check_delVD_bnds==1
    if transorbit.delVD*constants.normV2trueV > filters.delVD_max || transorbit.delVD*constants.normV2trueV < filters.delVD_min
        feas(9)=0;
    else
        feas(9)=1;
    end
else
    feas(9)=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check for N bounds
if filters.check_N_bnds==1
    if transorbit.N< filters.Nmin_bnd || transorbit.N> filters.Nmax_bnd
        feas(10)=0;
    else
        feas(10)=1;
    end
else
    feas(10)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if departure delv is less than 11km/s
% remember delVA is normalized, so bring it to original km/s
if filters.check_delVA_bnds==1
    if transorbit.delVA*constants.normV2trueV > filters.delVA_max || transorbit.delVA*constants.normV2trueV < filters.delVA_min
        feas(11)=0;
    else
        feas(11)=1;
    end
else
    feas(11)=1;
end
