function  F=optim_ulp(delv,Xrvmc,Wmc,rfstar,Tvec)
rfstar=rfstar(:);
F=0;
ERR=zeros(size(Xrvmc,1),3);
for i=1:size(Xrvmc,1)
    X=propagate_orbit_rv_mu([Xrvmc(i,1),Xrvmc(i,2),Xrvmc(i,3),Xrvmc(i,4)+delv(1),Xrvmc(i,5)+delv(2),Xrvmc(i,6)+delv(3)],Tvec,1);
   err=rfstar-X(end,1:3)';
   
%    ERR(i,:)=err;
%    C=cov(err);
%    F= F+ Wmc(i)*(err'*err);


    F= F+ Wmc(i)*(err'*err);
end

% F=det(cov(ERR));





