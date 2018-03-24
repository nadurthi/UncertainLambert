
function plot_earth_sphere(Xsphere,ToPo)
% function to rotate earth to time ET


% if ET is NaN : do not use earth rotation from J2000 frame
% These have to be loaded apriori
% ToPo=load('topo.mat','topo','topomap1');
% Xsphere.radii=cspice_bodvrd( 'EARTH', 'RADII', 3);
% SPICE kernels have to be loaded
% [x,y,z] = sphere(50);
% x=-x;
% y=-y;
from_frame='J2000';
to_frame='ITRF93';

props.AmbientStrength           = 0.1;
props.DiffuseStrength           = 1;
props.SpecularColorReflectance  = .5;
props.SpecularExponent          = 20;
props.SpecularStrength          = 1;
props.FaceColor                 = 'texture';
props.EdgeColor                 = 'none';
props.FaceLighting              = 'phong';
props.Cdata                     = ToPo.topo;
 
colormap(ToPo.topomap1);

x=zeros(size(Xsphere.x) );
y=zeros(size(Xsphere.y) );
z=zeros(size(Xsphere.z) );
R=eye(3);

for i=1:1:size(Xsphere.x,1)
    for j=1:1:size(Xsphere.x,2)        
        P= R*[Xsphere.x(i,j);Xsphere.y(i,j);Xsphere.z(i,j)];
        x(i,j)=P(1);y(i,j)=P(2);z(i,j)=P(3);
    end
end
radii=Xsphere.radii;
surface(x*radii(1),y*radii(2),z*radii(3),props);
% set(gca,'color','black')
alpha 0.5
% axis equal
grid on

hold on

% plot J2000 axis
% p=1.05*radii(1);
% plot3([0,p],[0,0],[0,0],'r--','linewidth',2)
% plot3([0,0],[0,p],[0,0],'y--','linewidth',2)
% plot3([0,0],[0,0],[0,p],'g--','linewidth',2)


% % plot ITRF93 coordinate axis
% X=[p,0,0;
%     0,p,0;
%     0,0,p];
% 
% 
% for i=1:3
%    X(i,:)=R*X(i,:)'; 
%    G=[0,0,0;X(i,:)];
%    plot3(G(:,1),G(:,2),G(:,3),'r--','linewidth',1)
%    plot3(G(:,1),G(:,2),G(:,3),'y--','linewidth',1)
%    plot3(G(:,1),G(:,2),G(:,3),'--','linewidth',1)
% end


end
