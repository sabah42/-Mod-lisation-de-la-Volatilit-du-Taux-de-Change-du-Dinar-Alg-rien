function [x] = MAR_ARCH( phi0,phi,alph,beta0,beta,n)
[p,K]=size(phi);
[q,K]=size(beta);
 x=zeros(2*n,1);
a0=sum(alph.*phi0);
a=1;
for i=1:p
    for k=1:K
        a=a-(phi(i,k)*alph(k));
    end
end
x(1:p,1)=a0/a;
for t=p+q+1:2*n
    k=randsample(K,1,true,alph);
    x(t)=phi0(k);
    for i=1:p
        x(t)=x(t)+phi(i,k)*x(t-i);
    end
% end
% for t=p+q+1:n
    h(t,k)=beta0(k);
    for j=1:q
        h(t,k)=h(t,k)+beta(j,k)*x(t-j)^2;
    end 
    x(t)=x(t)+sqrt(h(t,k))*randn;       
end
x=x(n+1:2*n,1);
end
