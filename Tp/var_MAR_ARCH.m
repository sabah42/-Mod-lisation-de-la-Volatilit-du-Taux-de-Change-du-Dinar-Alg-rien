function[var]=var_MAR_ARCH(Alpha,Phi0,Phi,Beta0,Beta,x,vectp,vectq)
[p,K]=size(Phi);
[q,K]=size(Beta);
n=length(x);
s1=zeros(n);
s2=zeros(n);
s3=zeros(n);
for t=p+1:n
    for k=1:K
        mu(t,k)=Phi0(k);
        for i=1:vectp(k)
            mu(t,k)=mu(t,k)+Phi(i,k)*x(t-i);
        end
    end
end
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
    end
end
for k=1:K
    for t=p+q+1:n
        s1(t)=s1(t)+Alpha(k)*h(t,k);
        s2(t)=s2(t)+Alpha(k)*mu(t,k)^2;
        s3(t)=s3(t)+Alpha(k)*mu(t,k);
    end
end
for t=p+q+1:n
    var(t)=s1(t)+s2(t)-(s3(t))^2;
end
end
