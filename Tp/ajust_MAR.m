function[x_ajust,eps_ajust,x]=ajust_MAR(x,Alpha,Phi0,Phi,Sigma,vect)
n=size(x,1);
[p,K]=size(Phi);
for t=p+1:n
    for k=1:K
        e(t,k)=x(t)-Phi0(k);
        for i=1:vect(k)
            e(t,k)=e(t,k)-Phi(i,k)*x(t-i);
        end
        s(t,k)=(Alpha(k)/Sigma(k))*(normpdf(e(t,k)/Sigma(k)));
    end
    for k=1:K
        tau(t,k)=(s(t,k)/sum(s(t,:)));
    end
end
for t=p+1:n
    y=max(tau(t,:));
    k=find(y==tau(t,:));
    x_ajust(t)=Phi0(k);
    for i=1:p
        x_ajust(t)=x_ajust(t)+Phi(i,k)*x(t-i);
        eps_ajust(t)=(x(t)-x_ajust(t))*inv(sqrt(Sigma(k)));
    end
end
end



