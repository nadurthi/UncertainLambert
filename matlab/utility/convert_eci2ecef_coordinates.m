function Y=convert_eci2ecef_coordinates(ET,X,constants)
from_frame='J2000';
to_frame='ITRF93';


Y=zeros(size(X(:,1:3)));
for i=1:length(ET)
    
    R=cspice_pxform( from_frame, to_frame, ET(i) );
    Y(i,:)=R*X(i,1:3)';
    
end

% [~,~,Alt]=cspice_recgeo(Y',constants.radii(1),constants.f);
