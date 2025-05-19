function [alp,tau,e,h,p,q,K] = E_step( phi0,phi,alph,beta0,beta,x,vectp,vectq)
n=size(x,1);
[p,K]=size(phi);
[q,K]=size(beta);
tau=zeros(n,K);
e=zeros(n,K);
h=zeros(n,K);
stau=zeros(1,K);
s=zeros(n,K);
%E_step
for t=p+q+1:n
    for k=1:K
        for i=1:vectp(k)
            e(t,k)=e(t,k)-phi(i,k)*x(t-i);
        end
        e(t,k)=e(t,k)+x(t)-phi0(k);
        for j=1:vectq(k)  
            h(t,k)=h(t,k)+beta(j,k)*e(t-j,k)^2;
        end
        h(t,k)=h(t,k)+beta0(k);
     s(t,k)=(alph(k)/sqrt(h(t,k)))*(normpdf(e(t,k)/sqrt(h(t,k))));
    end 
    for k=1:K
        tau(t,k)=(s(t,k)/sum(s(t,:)));
    end
end
%M_step
for k=1:K
    stau(k)=sum(tau(:,k));
    alp(k)=stau(k)/(n-p-q);
end

end