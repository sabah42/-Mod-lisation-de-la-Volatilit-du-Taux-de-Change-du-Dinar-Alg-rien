function [para,BIC]=critere_ARCH(Phi0,Phi,Beta0,Beta,Alpha,x,vectp,vectq)
[p,K]=size(Phi);
[q,K]=size(Beta);
n=length(x);
N=n-p-q;
%***************** le nbr de paramètres estimés non nul ******************* 
para=0;
for k=1:K
    if Alpha(k)~=0
        para=para+1;
    end
    if Beta0(k)~=0
        para=para+1;
    end
    for i=1:vectq(k)
        if Beta(i,k)~=0
            para=para+1;
        end
    end
   if Phi0(k)~=0
        para=para+1;
   end
    for i=1:vectp(k)
        if Phi(i,k)~=0
            para=para+1;
        end
    end
end
%************************ Le critère BIC **********************************
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
end
for t=p+q+1:n
for k=1:K
    l(t)=sum(s(t,:));
end
 L=sum(log(l(t)))/N;
 BIC=-2*N*L+log(n-p-q)*(para-1);
end
end