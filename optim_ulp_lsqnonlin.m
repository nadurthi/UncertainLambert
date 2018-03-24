function  F=optim_ulp_lsqnonlin(delv,Xrvmc,Wmc,rfstar,Tvec)
N=size(Xrvmc,1);
rfstar=rfstar(:);
F=zeros(N,1);
for i=1:N
    X=propagate_orbit_rv_mu([Xrvmc(i,1),Xrvmc(i,2),Xrvmc(i,3),Xrvmc(i,4)+delv(1),Xrvmc(i,5)+delv(2),Xrvmc(i,6)+delv(3)],Tvec,1);
    err=rfstar-X(end,1:3)';
    F(i) = Wmc(i)*sqrt((err'*err));
end






