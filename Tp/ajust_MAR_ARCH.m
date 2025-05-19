function[x_ajust,eps_ajust]=ajust_MAR_ARCH(x,Alpha,Phi0,Phi,Beta0,Beta,vectp,vectq)
n=size(x,1);
[p,K]=size(Phi);
[q,K]=size(Beta);
for t=p+1:n
    for k=1:K
        e(t,k)=x(t)-Phi0(k);
        for i=1:vectp(k)
            e(t,k)=e(t,k)-Phi(i,k)*x(t-i);
        end
    end
end

for t=p+q+1:n
    for k=1:K
        h(t,k)=Beta0(k);
        for j=1:vectq(k)
            h(t,k)=h(t,k)+Beta(j,k)*e(t-j,k)^2;
        end
        s(t,k)=(Alpha(k)/sqrt(h(t,k)))*(normpdf(e(t,k)/sqrt(h(t,k))));
    end
 for k=1:K   
    tau(t,k)=(s(t,k)/sum(s(t,:)));
 end
end
for t=p+q+1:n
    y=max(tau(t,:));
    k=find(y==tau(t,:));
    x_ajust(t)=Phi0(k);
    for i=1:p
        x_ajust(t)=x_ajust(t)+Phi(i,k)*x(t-i);
        eps_ajust(t)=(x(t)-x_ajust(t))*inv(sqrt(h(t,k)));
    end
end
end
