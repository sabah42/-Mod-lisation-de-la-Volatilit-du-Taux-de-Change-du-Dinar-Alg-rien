function [ddphi] = ddphiest(dedphi,dhdphi,tau,h,p,q,vectp,vectq)
[n,K]=size(tau);
ddphi=zeros(p+1,p+1,K);
for k=1:K
for t=p+q+1:n
ddphi(1,1,k)=ddphi(1,1,k)+(tau(t,k)*(dhdphi(1,t,k)^2/(2*h(t,k)^2))+(dedphi(1,t,k)^2/h(t,k)));
for i=1:vectp(k)
   ddphi(i+1,1,k)=ddphi(i+1,1,k)+(tau(t,k)*(dhdphi(i+1,t,k)*dhdphi(1,t,k)/(2*h(t,k)^2))+(dedphi(1,t,k)*dedphi(i+1,t,k)/h(t,k))); 
   for j=1:vectq(k)
ddphi(i+1,j+1,k)=ddphi(i+1,j+1,k)+(tau(t,k)*(dhdphi(i+1,t,k)*dhdphi(j+1,t,k)/(2*h(t,k)^2))+(dedphi(j+1,t,k)*dedphi(i+1,t,k)/h(t,k)));
   end
end
for j=1:vectq(k)
    ddphi(1,j+1,k)=ddphi(1,j+1,k)+(tau(t,k)*(dhdphi(1,t,k)*dhdphi(j+1,t,k)/(2*h(t,k)^2))+(dedphi(j+1,t,k)*dedphi(1,t,k)/h(t,k)));
end

for i=1:vectp(k)
    for j=1:vectq(k)
        ddphi(i+1,j+1)=-(ddphi(i+1,j+1)) *inv((n-p-q));
    end
end
end
end
end