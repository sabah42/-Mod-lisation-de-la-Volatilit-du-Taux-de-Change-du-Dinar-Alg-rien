function [dphi] = dphiest( dedphi,dhdphi,tau,e,h,p,q,vectp)
[n,K]=size(tau);
dphi=zeros(p+1,K);
for k=1:K
     for t=p+q+1:n
           dphi(1,k)=dphi(1,k)+((tau(t,k)/(2*h(t,k)))*dhdphi(1,t,k)*((e(t,k)^2/h(t,k))-1)-(tau(t,k)*e(t,k)*dedphi(1,t,k)/(h(t,k))))*inv(n-p-q);
     end
    
    for i=1:vectp(k)
        for t=p+q+1:n
           dphi(i+1,k)=dphi(i+1,k)+((tau(t,k)/(2*h(t,k)))*dhdphi(i+1,t,k)*((e(t,k)^2/h(t,k))-1)-(tau(t,k)*e(t,k)*dedphi(i+1,t,k)/(h(t,k))))*inv(n-p-q);
        end
    end
end
 
end
