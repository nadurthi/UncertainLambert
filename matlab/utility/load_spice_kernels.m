function load_spice_kernels(p)

comp=computer();

if strcmp(comp,'GLNXA64')
   pth=[p,'/MATLAB/LINUX64_mice']; 
elseif  strcmp(comp,'PCWIN64')
    pth=[p,'/MATLAB/WINDOWS64_mice'];
elseif  strcmp(comp,'MACI64')
    pth=[p,'/MATLAB/MAC64_mice'];
end
addpath( [pth,'/lib']) % location of lib folder
addpath( [pth,'/src/mice']) 
addpath( [p,'/kernels']) % folder of required kernels

%  leap seconds kernel for Epephermis time
cspice_furnsh( [p,'/kernels/naif0012.tls'] )

% kernel for earth fixed frames
cspice_furnsh( [p,'/kernels/earth_070425_370426_predict.bpc'] )

% kernel for earth properties like radii
cspice_furnsh( [p,'/kernels/pck00010.tpc'] )



