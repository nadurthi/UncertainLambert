function [R,V]=normalize_rv(r,v,constants)

R=r(:,1:3)./constants.RU ;
V=v(:,1:3).*(constants.TU/constants.RU); 