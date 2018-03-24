function update_handles=plot_ground_track(update_handles,ET,Xsphere,SatPos,SatTraj,SatTraj_ET,ToPo,plotprops)
% SatPos: position of satellite in J2000 frame [x;y;z]
% SatTraj can be NaN if you do not want to plot the trajectory. SatTraj_ET
%       is the ET time of the trajectory 
% plotprops.sat_marker='ro'
% plotprops.sattraj_marker='r--'
% plotprops.sattraj_line='r--';
% plotprops.MarkerFaceColor='r';

% update_handles : the plot handles of satellite position and trajectory.
%                   this is used for fast plotting and avoiding plottin the
%                   earth topography repeatedly.

from_frame='J2000';
to_frame='ITRF93';
R=cspice_pxform( from_frame, to_frame, ET );



SatPos=SatPos(:);


f = (Xsphere.radii(1)-Xsphere.radii(3))/Xsphere.radii(1);

SatPos_ITRF=R*SatPos;
[Satlong,Satlat,Satalt]=cspice_recgeo(SatPos_ITRF,Xsphere.radii(1),f);
if Satlong<=0
    Satlong=Satlong+2*pi;
end

if isnan(SatTraj)==0
    for i=1:size(SatTraj,1)
        R=cspice_pxform( from_frame, to_frame, SatTraj_ET(i) );
        SatTraj_ITRF=R*SatTraj(i,:)';
        [ln,lt,alt]=cspice_recgeo(SatTraj_ITRF,Xsphere.radii(1),f);
        if ln<=0
            ln=ln+2*pi;
        end
        SatTraj(i,:)=[ln,lt,alt];
    end
    
    
end

% image([0 360],[-90 90], ToPo.topo(end:-1:1,end:-1:1), 'CDataMapping', 'scaled')
% colormap(ToPo.topomap1)
% hold on
% grid on

% [Satlat*180/pi,Satlong*180/pi]

for i=1:length(update_handles)

        delete(update_handles{i});

end
% keyboard
update_handles={};
update_handles{1}=plot(Satlong*180/pi,Satlat*180/pi,plotprops.sat_marker,'MarkerSize',plotprops.sat_markersize,'linewidth',2,'MarkerFaceColor',plotprops.MarkerFaceColor,'tag','satpos');
k=1;

if isnan(SatTraj)==0
    % first break the trajectories that cross from jump from 2*pi long to 0
    % or vice versa
    % and pi/2 to -pi/2 and vice versa
    meangap=mean( sqrt(sum(diff(SatTraj(:,[1,2])).^2,2)) ) ;
    k=1;
    p=1;
    SatTraj_parts={};
    for i=2:size(SatTraj,1)
        if  norm(SatTraj(i,[1,2])-SatTraj(i-1,[1,2])) >3*meangap
            SatTraj_parts{p}=SatTraj(k:i-1,[1,2]);
            p=p+1;
            k=i;
            
        end
    end
    
    SatTraj_parts{p}=SatTraj(k:end,[1,2]);
    
    
    for i= 1:length(SatTraj_parts)
        update_handles{k+i}=plot(SatTraj_parts{i}(:,1)*180/pi,SatTraj_parts{i}(:,2)*180/pi,plotprops.sattraj_line,'linewidth',2);
    end
end

axis equal                                % set axis units to be the same size

ax = gca;                                 % get current axis
ax.XLim = [0 360];                        % set x limits
ax.YLim = [-90 90];                       % set y limits
ax.XTick = [0 60 120 180 240 300 360];    % define x ticks
ax.YTick = [-90 -60 -30 0 30 60 90];      % define y ticks