function plot_transfer_static(fig,ax,perturbed_efmoutput,depind,tfind, minsol,propagator,constants,config)
% @brief: plot the animation of a given transfer orbit.
% @Input: fig - figure object
%         ax - axis object
%         depind - the departure index of the cell SimSol_transferorbits (the output from efmoutput)
%         tfind - the transfer_time index of the cell SimSol_transferorbits (the output from efmoutput)
%         minsol - the solution index for the cell {depind,tfind}
%         efmoutput - efmoutput the solution output of EFMcompute_kepler
% @Output: No output return, the figure is automatically updated at fixed
%           intervals

inputconfig=perturbed_efmoutput.inputconfig;

transferorbits=perturbed_efmoutput.SimSol_transferorbits_kepler{depind,tfind};
pert_transferorbits=perturbed_efmoutput.SimSol_transferorbits{depind,tfind};

transferorbit=convert_transferorbit_mat2struct(transferorbits(minsol,:));
pert_transferorbit=convert_transferorbit_mat2struct_perturbed(pert_transferorbits(minsol,:));
    
% transferorbit=efmoutput.SimSol_transferorbits{depind,tfind}{minsol};
% transferorbit;
% Now compute trajectories for the departure, transfer and arrival
dep_epoch=inputconfig.dep_epoch;
arr_epoch=inputconfig.arr_epoch;



dep_time=inputconfig.DepTime(depind) ;
transfer_time=inputconfig.TimeFlights(tfind);

starttime=0.1*max([dep_epoch,arr_epoch])+0.9*dep_time ;

endtime=dep_time+transfer_time;

simtime_to_deptime=linspace(starttime,dep_time,1000);
simtime_to_arrtime=linspace(dep_time,endtime,1000);


% simulate to departure time
Tvec=[dep_epoch,linspace(dep_epoch,dep_epoch+inputconfig.DepP,3000)];
dep_orbit_1period=propagate_orbit_rv([inputconfig.rD,inputconfig.vD],Tvec,constants);

Tvec=[arr_epoch,linspace(arr_epoch,arr_epoch+inputconfig.ArrP,3000)];
arr_orbit_1period=propagate_orbit_rv([inputconfig.rA,inputconfig.vA],Tvec,constants);


Tvec=[dep_epoch,simtime_to_deptime];
dep_orbit=propagate_orbit_rv([inputconfig.rD,inputconfig.vD],Tvec,constants);
dep_orbit_phase1=dep_orbit(2:end,:);

Tvec=[arr_epoch,simtime_to_deptime];
arr_orbit=propagate_orbit_rv([inputconfig.rA,inputconfig.vA],Tvec,constants);
arr_orbit_phase1=arr_orbit(2:end,:);

Tvec=[simtime_to_deptime(end),simtime_to_arrtime];
arr_orbit=propagate_orbit_rv(arr_orbit_phase1(end,:),Tvec,constants);
arr_orbit_phase2=arr_orbit(2:end,:);

Tvec=[simtime_to_deptime(end),simtime_to_arrtime]-simtime_to_deptime(end);
tran_orbit=propagate_orbit_rv([transferorbit.r1*constants.normX2trueX,transferorbit.v1*constants.normV2trueV],Tvec,constants);
tran_orbit_phase=tran_orbit(2:end,:);


Tvec=[simtime_to_arrtime];
[t,pert_tran_orbit]=propagator(Tvec,[pert_transferorbit.r1*constants.normX2trueX,pert_transferorbit.v1*constants.normV2trueV]);
pert_tran_orbit_phase=pert_tran_orbit;




plotprops.sat_marker='ro';
plotprops.sat_markersize=7;
plotprops.sattraj_line='r--';
plotprops.MarkerFaceColor='r';


% 
% for i=1:length(simtime_to_deptime)
%     figure(fig)
% 
%     if config.useSPICE
%         subplot(1,2,1)
%     end
%     hold off;
%     plot3(dep_orbit_1period(:,1),dep_orbit_1period(:,2),dep_orbit_1period(:,3),'r')
%     grid on
%     hold on
%     plot3(arr_orbit_1period(:,1),arr_orbit_1period(:,2),arr_orbit_1period(:,3),'b')
%     
%     
%     plot3(dep_orbit_phase1(i,1),dep_orbit_phase1(i,2),dep_orbit_phase1(i,3),'ro','linewidth',2,'MarkerSize',6)
%     plot3(arr_orbit_phase1(i,1),arr_orbit_phase1(i,2),arr_orbit_phase1(i,3),'bs','linewidth',2,'MarkerSize',6)
%     
%     
% 
%     if config.useSPICE
%         subplot(1,2,2)
%         update_handles=plot_ground_track(update_handles,simtime_to_deptime(i),Xsphere,dep_orbit_phase1(i,1:3),NaN,NaN,ToPo,plotprops);
%     end
%     
%     pause(0.001)
% %     saveas(gcf,sprintf('%05d',i),'jpg')
% end
% k=i;
% update_handles={};

    figure(fig)

    plot3(dep_orbit_1period(:,1),dep_orbit_1period(:,2),dep_orbit_1period(:,3),'r')
    grid on
    hold on
    plot3(arr_orbit_1period(:,1),arr_orbit_1period(:,2),arr_orbit_1period(:,3),'b')
    
    
    plot3(tran_orbit_phase(:,1),tran_orbit_phase(:,2),tran_orbit_phase(:,3),'k','linewidth',2,'MarkerSize',6)
    plot3(pert_tran_orbit_phase(:,1),pert_tran_orbit_phase(:,2),pert_tran_orbit_phase(:,3),'g--','linewidth',2,'MarkerSize',6)

%     plot3(pert_tran_orbit_phase(:,1),pert_tran_orbit_phase(:,2),pert_tran_orbit_phase(:,3),'g','linewidth',2,'MarkerSize',6)
    
%     plot3(arr_orbit_phase2(:,1),arr_orbit_phase2(:,2),arr_orbit_phase2(:,3),'b','linewidth',2,'MarkerSize',6)
    
    
    
legend('Departure orbit','Arrival orbit','Transfer Orbit','Perturbed Transfer Orbit')
