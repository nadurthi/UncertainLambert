function [v,ep]=linear_robust_lambert(PHIs,Vs)

I3=eye(3);
A=[];
B=[];
for i=1:1:length(PHIs)
    phi12i=PHIs{i}(1:3,4:6);
    vi=Vs{i};
    vi=vi(:);
    a=[ phi12i,-I3;-phi12i,-I3];
    b=[ phi12i*vi; -phi12i*vi];
    
    A=vertcat(A,a);
    B=vertcat(B,b);
    
end

f=[0;0;0;1;1;1];
x = linprog(f,A,B,[],[],[-inf,-inf,-inf,0,0,0]);

v=x(1:3);
ep=x(4:6);








