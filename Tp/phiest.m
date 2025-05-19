function [Phi,ph0] = phiest( phi0,phi,beta,tau,e,h,x,p,q,K,vectp,vectq)
n=size(x,1);
dphi=zeros(p+1,K);
ddphi=zeros(p+1,p+1,K);
phii=[phi0;phi];
for k=1:K
    for i=1:vectp(k)+1
        if i==1
            for t=p+q+1:n
                l=zeros(n,K);
                for z=1:vectq(k)
                    l(t,k)=l(t,k)+beta(z,k)*e(t-z,k);
                end
                dphi(i,k)=dphi(i,k)+((tau(t,k)/h(t,k))*(-l(t,k))*(((e(t,k)^2)/(h(t,k)))-1)+(tau(t,k)*(e(t,k))/(h(t,k))))/(n-p-q);
            end
            for j=1:vectp(k)+1
                if j==1
                    for t=p+q+1:n
                        l=zeros(n,K);
                        for z=1:vectq(k)
                            l(t,k)=l(t,k)+beta(z,k)*e(t-z,k);
                        end
                        ddphi(i,j,k)=ddphi(i,j,k)+(tau(t,k)*((2*l(t,k)^2/h(t,k)^2)+1/h(t,k)))/(n-p-q);
                    end
                else
                    
                    for t=p+q+1:n
                        l=zeros(n,K);
                        for z=1:vectq(k)
                            l(t,k)=l(t,k)+beta(z,k)*e(t-z,k);
                        end
                        r=zeros(n,K);
                        for z=1:vectq(k)
                            r(t,k)=r(t,k)+beta(z,k)*e(t-z,k)*x(t-j-z+1);
                        end
                        ddphi(i,j,k)=ddphi(i,j,k)+(tau(t,k)*(2*l(t,k)*r(t,k))/(h(t,k)^2)+x(t-j+1)/h(t,k))/(n-p-q);
                    end
                end
            end
        else
            for t=p+q+1:n
                r=zeros(n,K);
                for z=1:vectq(k)
                    r(t,k)=r(t,k)+beta(z,k)*e(t-z,k)*x(t-z-i+1);
                end
                dphi(i,k)=dphi(i,k)+((tau(t,k)/h(t,k))*(-r(t,k))*((e(t,k)^2/h(t,k))-1)+(tau(t,k)*e(t,k)*x(t-i+1)/h(t,k)))/(n-p-q);
            end
            for j=1:vectp(k)+1
                if j==1
                    for t=p+q+1:n
                        l=zeros(n,K);
                        for z=1:vectq(k)
                            l(t,k)=l(t,k)+beta(z,k)*e(t-z,k)*x(t-i-z+1);
                        end
                        r=zeros(n,K);
                        for z=1:vectq(k)
                            r(t,k)=r(t,k)+beta(z,k)*e(t-z,k);
                        end
                        ddphi(i,j,k)=ddphi(i,j,k)+(tau(t,k)*(2*l(t,k)*r(t,k))/(h(t,k)^2)+x(t-i+1)/h(t,k))/(n-p-q);
                    end
                else
                    for t=p+q+1:n
                        l=zeros(n,K);
                        for z=1:vectq(k)
                            l(t,k)=l(t,k)+beta(z,k)*e(t-z,k)*x(t-i-z+1);
                        end
                        r=zeros(n,K);
                        for z=1:vectq(k)
                            r(t,k)=r(t,k)+beta(z,k)*e(t-z,k)*x(t-j-z+1);
                        end
                        ddphi(i,j,k)=ddphi(i,j,k)+(tau(t,k)*(2*l(t,k)*r(t,k))/(h(t,k)^2)+x(t-j+1)*x(t-i+1)/h(t,k))/(n-p-q);
                    end
                end
            end
        end
    end
   
   ddphi(:,:,k)=inv(ddphi(:,:,k));
end

for k=1:K
    for i=1:vectp(k)+1
        for j=1:vectp(k)+1
          
            PHII(i,k)=phii(i,k)+(dphi(i,k)*ddphi(i,j,k));
           
        end
        
    end
    
end
 ph0=PHII(1,:);
for k=1:K
    for i=1:vectp(k)
        Phi(i,k)=PHII(i+1,k);
    end
end
end







