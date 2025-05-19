function [f] = fct_predic( phi0,phi,sigma,alph,t,x,y)
n=size(x,1);
[p,K]=size(phi);
f=0;
for k=1:K
    e(k)=y-phi0(k);
    for i=1:p
        e(k)=e(k)-phi(i,k)*x(t-i);
    end
    
   f=f+alph(k)*inv(sigma(k))*normpdf(e(k)*inv(sigma(k)));
end
end
