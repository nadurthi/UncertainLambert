function [TimeFlights_normalized,DepOrbit_normalized,ArrOrbit_normalized_cell]=GetDepArr_normalized_orbit_positions_perturbed(inputconfig,propagator,constants,config)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using the grid of departure times: compute the corresponding grid of
% positions along the departure orbit by propagating the epoch element set/
% coordinates. Now using the depature time grid and Transfer time grid,
% compute the corresponding arrival positions on the arrival orbit
% Finally, normalize all data using constants
%       normalized time     -> true time/constants.TU
%       normalized position -> true position/constants.RU
%       normalized velocity -> true velocity*(constants.TU/constants.RU



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% departure orbit
R    = inputconfig.rD(:);
V    = inputconfig.vD(:);

Tvec = [inputconfig.dep_epoch,inputconfig.DepTime];
[t,DepOrbit]=propagator(Tvec,[R;V]);
% keyboard
[DepOrbit_normalized(:,1:3),DepOrbit_normalized(:,4:6)]=normalize_rv(DepOrbit(2:end,1:3),DepOrbit(2:end,4:6),constants);

%% for all departure time and transfer time, compute the arrival state

R    = inputconfig.rA;
V    = inputconfig.vA;
ArrOrbit_normalized_cell=cell(length( inputconfig.DepTime ));


parfor i=1:1:length( inputconfig.DepTime )
    ArrOrbit_normalized_cell{i}=cell(1,length(inputconfig.TimeFlights));
%     S1=cell(1,length(inputconfig.TimeFlights));
    for j=1:1:length(inputconfig.TimeFlights)
        [i,j]
        T=inputconfig.DepTime(i)+inputconfig.TimeFlights(j);
        Tvec = [inputconfig.arr_epoch,T];
        [t,ArrOrbit]=propagator(Tvec,[R,V]);
        [r,v]=normalize_rv(ArrOrbit(2:end,1:3),ArrOrbit(2:end,4:6),constants);
        ArrOrbit_normalized_cell{i}{j}=[r(end,:),v(end,:)];
    end
%         for j=1:length(inputconfig.TimeFlights)
%             =S1{j};
%         end
end


TimeFlights_normalized     = inputconfig.TimeFlights'./constants.TU;

