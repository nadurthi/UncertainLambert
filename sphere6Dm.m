function X1 = sphere6Dm(n)

N=n^5;
X1 = zeros(N,n);

% ph1,ph2,ph3,ph4,ph5
ph1 =linspace(10*pi/180,170*pi/180,n);
ph2 =linspace(10*pi/180,170*pi/180,n);
ph3 =linspace(10*pi/180,170*pi/180,n);
ph4 =linspace(10*pi/180,170*pi/180,n);
ph5 =linspace(10*pi/180,350*pi/180,n);

a=1;
b=1;

for i=1:n
    for j=1:n
        for k=1:n
            for p=1:n
                for q=1:n
                    X1(a,1)=cos(ph1(i));
                    X1(a,2)=sin(ph1(i))*cos(ph2(j));
                    X1(a,3)=sin(ph1(i))*sin(ph2(j))*cos(ph3(k));
                    X1(a,4)=sin(ph1(i))*sin(ph2(j))*sin(ph3(k))*cos(ph4(p));
                    X1(a,5)=sin(ph1(i))*sin(ph2(j))*sin(ph3(k))*sin(ph4(p))*cos(ph5(q));
                    X1(a,6)=sin(ph1(i))*sin(ph2(j))*sin(ph3(k))*sin(ph4(p))*sin(ph5(q));
%                     X1(a,:)
%                     if norm(X1(a,:))<0.999 || norm(X1(a,:))>1.001
%                         X1(a,:)
%                     end
                    a=a+1;
                end
            end
        end
    end
    
%     keyboard
    
    
end


% H=10^6*X1(:,1)+10^5*X1(:,2)+10^4*X1(:,3)+10^3*X1(:,5)+10^2*X1(:,5)+X1(:,6);
% Hunq=unique(H);

X1=unique(X1,'rows');








end