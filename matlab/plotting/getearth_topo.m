function [Xsphere,ToPo]=getearth_topo(radii)

ToPo=load('topo.mat','topo','topomap1');
Xsphere.radii=radii;
[x,y,z] = sphere(50);
x=-x;
y=-y;

Xsphere.x=x;
Xsphere.y=y;
Xsphere.z=z;
