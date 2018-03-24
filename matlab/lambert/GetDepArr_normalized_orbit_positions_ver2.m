function [TimeFlights_normalized,DepOrbit_normalized,ArrOrbit_normalized_cell]=GetDepArr_normalized_orbit_positions_ver2(inputconfig,constants,mode)
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
R    = inputconfig.rD;
V    = inputconfig.vD;

Tvec = [inputconfig.dep_epoch,inputconfig.DepTime]-inputconfig.dep_epoch;
DepOrbit=propagate_orbit_rv([R,V],Tvec,constants);
[DepOrbit_normalized(:,1:3),DepOrbit_normalized(:,4:6)]=normalize_rv(DepOrbit(2:end,1:3),DepOrbit(2:end,4:6),constants);

% rr=[0.0919088536398892653;
%   -1.04473493915354276;
% 0.00144205432889136262];
% vv=[0.603720476521701799;
%  0.0546169279488273457;
%   0.765815924582433416];
% [ r, v, Ehat ] = FnG(0, 100, R*constants.trueX2normX, V*constants.trueV2normV, 1)

%% for all departure time and transfer time, compute the arrival state

R    = inputconfig.rA;
V    = inputconfig.vA;
ArrOrbit_normalized_cell=cell(length( inputconfig.DepTime ),length( inputconfig.TimeFlights ));


for i=1:1:length( inputconfig.DepTime )
    for j=1:1:length(inputconfig.TimeFlights)
        T=inputconfig.DepTime(i)+inputconfig.TimeFlights(j);
            Tvec = [inputconfig.arr_epoch,T]-inputconfig.arr_epoch;
            ArrOrbit=propagate_orbit_rv([R,V],Tvec,constants);
            [r,v]=normalize_rv(ArrOrbit(2:end,1:3),ArrOrbit(2:end,4:6),constants);
            ArrOrbit_normalized_cell{i,j}={[r(end,:),v(end,:)],T};
    end
end


TimeFlights_normalized     = inputconfig.TimeFlights'./constants.TU;

