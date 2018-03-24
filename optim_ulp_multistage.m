function  F=optim_ulp_multistage(alphas,Xrvmc,Wmc,rfstar,Tvec)
rfstar=rfstar(:);
F=0;
for i=1:size(Xrvmc,1)
    X=Xrvmc(i,:);
    for ti=1:1:length(Tvec)-1
        X=propagate_orbit_rv_mu([X(1),X(2),X(3),X(4)+alphas(3*(ti-1)+1),X(5)+alphas(3*(ti-1)+2),X(6)+alphas(3*(ti-1)+3)],[Tvec(ti),Tvec(ti+1)],1);
        X=X(end,:);
    end
    X=X(:)';
   err=rfstar-X(end,1:3)';
   F= F+ Wmc(i)*(err'*err);
end






