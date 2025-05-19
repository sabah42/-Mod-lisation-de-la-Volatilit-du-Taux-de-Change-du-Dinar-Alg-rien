function [alp,ph0,Ph,bet0,bet] = EM_step( phi0,phi,alph,beta0,beta,x,vectp,vectq)
n=length(x);
[p,K]=size(phi);
[q,K]=size(beta);
N=n-p-q;

%*********************** E_step *******************************************

for t=p+1:n
    for k=1:K
        e(t,k)=x(t)-phi0(k);
        for i=1:vectp(k)
            e(t,k)=e(t,k)-phi(i,k)*x(t-i);
        end
    end
end

for t=p+q+1:n
    for k=1:K
        h(t,k)=beta0(k);
        for j=1:vectq(k)
            h(t,k)=h(t,k)+beta(j,k)*e(t-j,k)^2;
        end
        s(k)=(alph(k)/sqrt(h(t,k)))*(normpdf(e(t,k)/sqrt(h(t,k))));
    end
    
    tau(t,:)=(s/sum(s));
end

%********************** M_step ********************************************

%********************** estimation de alpha *******************************

for k=1:K
    stau(k)=sum(tau(:,k));
    alp(k)=stau(k)/N;
end
% les dérivées
for k=1:K
    for t=p+1:n
        dedphi(1,t,k)=-1;
        for i=1:vectp(k)
            dedphi(i+1,t,k)=-x(t-i);
        end
    end
end
for k=1:K
    for t=p+q+1:n
        for j=1:vectq(k)
            E1(j)=beta(j,k)*e(t-j,k);
        end
        dhdphi(1,t,k)=-2*sum(E1);
        for i=1:vectp(k)
            for j=1:vectq(k)
                E2(j)=beta(j,k)*e(t-j,k)*dedphi(i+1,t-j,k);
            end
            dhdphi(i+1,t,k)=2*sum(E2);
        end
    end
end

%*********************** estimation de phi ********************************

for k=1:K
    phii(1,k)=phi0(k);
    for i=1:vectp(k)
        phii(i+1,k)=phi(i,k);
    end
end
[dphi] = dphiest( dedphi,dhdphi,tau,e,h,p,q,vectp);
[ddphi] = ddphiest(dedphi,dhdphi,tau,h,p,q,vectp,vectq);
for k=1:K
    PHII(:,k)=phii(:,k)-(inv(ddphi(:,:,k))*dphi(:,k));
    ph0(k)=PHII(1,k);
    for i=1:vectp(k)
        Ph(i,k)=PHII(i+1,k);
    end
end

%*********************** estimation de beta *******************************

for k=1:K
    betaa(1,k)=beta0(k);
    for i=1:vectq(k)
        betaa(i+1,k)=beta(i,k);
    end
end
[dbeta,ddbeta] = dbetest( ph0,Ph,beta,beta0,tau,x,vectp,vectq);
for k=1:K
    BETA(:,k)=betaa(:,k)-(inv(ddbeta(:,:,k))*dbeta(:,k));
    if BETA(1,1)>0
        bet0(k)=BETA(1,k);
    else
        bet0(k)=-BETA(1,k);
    end
    for i=1:vectq(k)
        if BETA(i+1,k)>=0
            bet(i,k)=BETA(i+1,k);
        else
            bet(i,k)=-BETA(i+1,k);
        end
    end
end

end

