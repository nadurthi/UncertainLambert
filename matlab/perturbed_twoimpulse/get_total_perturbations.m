function xyz_acc_pert = get_total_perturbations(t,x,perturb_models,constants,Satmodel,config)
% perturb_models is a cell of structs with the function call, the arguments
% required
xyz_acc_pert=[0;0;0];
for i=1:length(perturb_models)
    pertb_model=perturb_models{i};
    args={t,x,constants,Satmodel,config}; % load default arguments for all perturbation models
                                   % the default uses eci coordinates
    xyz_acc_pert=xyz_acc_pert+pertb_model(args{:});

end

