function [Alpha,Phi0,Phi,Beta0,Beta]=EM_MARARCH(x,alpha,phi0,phi,beta0,beta,itmax)
n=length(x);
K=size(alpha,2);
p=size(phi,1);
q=size(beta,1);
N=n-p-q;
for it=1:itmax
    
    for t=p+1:n
        for k=1:K
            e(t,k)=x(t)-phi0(it,k);
            for i=1:p
                e(t,k)=e(t,k)-phi(i,k,it)*x(t-i);
            end
        end
    end
    for t=p+q+1:n
        for k=1:K
            h(t,k)=beta0(it,k);
            for j=1:q
                h(t,k)=h(t,k)+beta(j,k,it)*((e(t-j,k))^2);
            end
            S(k)=alpha(it,k)*inv(sqrt(h(t,k)))*normpdf(e(t,k)*inv(sqrt(h(t,k))));
        end
        tau(t,:)=S/sum(S);
    end
    
    alpha(it+1,:)=sum(tau)/N;
    for k=1:K
        PHI(1,k)=phi0(it,k);
        for i=1:p
            PHI(i+1,k)=phi(i,k,it);
        end
    end
    for k=1:K
        for t=p+1:n
            dedphi(1,t,k)=-1;
            for i=1:p
                dedphi(i+1,t,k)=-x(t-i);
            end
        end
    end
    for k=1:K
        for t=p+q+1:n
            for j=1:q
                E1(j)=beta(j,k,it)*e(t-j,k);
            end
            dhdphi(1,t,k)=-2*sum(E1);
            for i=1:p
                for j=1:q
                    E2(j)=beta(j,k,it)*e(t-j,k)*dedphi(i+1,t-j,k);
                end
                dhdphi(i+1,t,k)=2*sum(E2);
            end
        end
    end
    for k=1:K
        
        for i=0:p
            for t=p+q+1:n
                S3(t)=tau(t,k)*(inv(2*h(t,k)))*dhdphi(i+1,t,k)*((e(t,k)^2*inv(h(t,k)))-1);
                S4(t)=tau(t,k)*e(t,k)*inv(h(t,k))*dedphi(i+1,t,k);
            end
            dldphi(i+1,k)=(1/N)*(sum(S3)-sum(S4));
        end
    end
    for k=1:K
        for i=0:p
            for j=0:p
                for t=p+q+1:n
                    A3(t)=tau(t,k)*(inv(2*h(t,k)^2))*dhdphi(i+1,t,k)*dhdphi(j+1,t,k)+(inv(h(t,k))*dedphi(i+1,t,k)*dedphi(j+1,t,k));
                end
                d2ldphiij(i+1,j+1,k)=-(1/N)*sum(A3);
            end
        end
    end
    for k=1:K
        PP(:,k)=PHI(:,k)-inv(d2ldphiij(:,:,k))*dldphi(:,k);
        phi0(it+1,k)=PP(1,k);
        for i=1:p
            phi(i,k,it+1)=PP(i+1,k);
        end
    end
    %estimation des beta%%%%%
    for k=1:K
        BETA(1,k)=beta0(it,k);
        for i=1:q
            BETA(i+1,k)=beta(i,k,it);
        end
    end
    %calcul de e et h
    for t=p+q+1:n
        for k=1:K
            for j=0:q
                e(t-j,k)=x(t-j)-phi0(it+1,k);
                for i=1:p
                    e(t-j,k)=e(t-j,k)-phi(i,k,it+1)*x(t-i-j);
                end
            end
            h(t,k)=beta0(it,k);
            for j=1:q
                h(t,k)=h(t,k)+beta(j,k,it)*(e(t-j,k)^2);
            end
        end
    end
    % les dérivées
    for k=1:K
        for t=p+q+1:n
            dhdbeta(1,t,k)=1;
            for i=1:q
                dhdbeta(i+1,t,k)=e(t-i,k)^2;
            end
        end
    end
    %estimation
    for k=1:K
        for i=0:q
            for t=p+q+1:n
                S6(t)=tau(t,k)*inv(2*h(t,k))*dhdbeta(i+1,t,k)*(((e(t,k)^2)*inv(h(t,k)))-1);
            end
            dldbeta(i+1,k)=(1/N)*sum(S6);
        end
    end
    for k=1:K
        for i=0:q
            for j=0:q
                for t=p+q+1:n
                    A6(t)=tau(t,k)*inv(2*(h(t,k)^2))*dhdbeta(i+1,t,k)*dhdbeta(j+1,t,k);
                end
                d2ldbetaij(i+1,j+1,k)=-(1/N)*sum(A6);
            end
        end
    end
    for k=1:K
        BB(:,k)=BETA(:,k)-(inv(d2ldbetaij(:,:,k)))*dldbeta(:,k);
        if BB(1,k)>0
            beta0(it+1,k)=BB(1,k);
        else
            beta0(it+1,k)=-BB(1,k);
        end
        if BB(i+1,k)>=0
            for i=1:q
                beta(i,k,it+1)=BB(i+1,k);
            end
        else
            for i=1:q
                beta(i,k,it+1)=-BB(i+1,k);
            end
        end
    end
end
    Alpha=alpha(end,:);
    Phi0=phi0(end,:);
    Phi=phi(:,:,end);
    Beta0=beta0(end,:);
    Beta=beta(:,:,end);
end
end










