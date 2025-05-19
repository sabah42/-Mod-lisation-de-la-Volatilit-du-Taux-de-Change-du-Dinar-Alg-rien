function [bet0,bet] = betest( ph0,Phi,beta,beta0,tau,x,vectq,vectp)
[n,K]=size(tau);
p=size(Phi,1);
q=size(beta,1);
dbeta=zeros(q+1,K);
ddbeta=zeros(q+1,q+1,K);
e=zeros(n,K);
h=zeros(n,K);
betaa=[beta0;beta];

for k=1:K
    for t=p+1:n
        e(t,k)=x(t)-ph0(k);
        for i=1:vectp(k)
            e(t,k)=e(t,k)-(Phi(i,k)*x(t-i));
        end
    end
    for t=p+q+1:n
        h(t,k)=beta0(k);
        for j=1:vectq(k)
            h(t,k)=h(t,k)+beta(j,k)*(e(t-j,k))^2;
        end
    end
end

for k=1:K
    for i=1:vectq(k)+1
        if i==1
            for t=p+q+1:n
                dbeta(i,k)=dbeta(i,k)+(tau(t,k)/(2*h(t,k)))*((e(t,k)^2/h(t,k))-1)/(n-p-q);
            end
            for j=1:vectq(k)+1
                
                for t=p+q+1:n
                    if j==1
                        ddbeta(i,j,k)=ddbeta(i,j,k)+(tau(t,k)/(2*h(t,k)^2))/(n-p-q);
                    else
                        ddbeta(i,j,k)=ddbeta(i,j,k)+((tau(t,k)/(2*h(t,k)^2))*e(t-j+1,k)^2)/(n-p-q);
                    end
                end
            end
        else
            for t=p+q+1:n
                dbeta(i,k)=dbeta(i,k)+((tau(t,k)/(2*h(t,k)))*((e(t,k)^2/h(t,k))-1)*e(t-i+1,k)^2)/(n-p-q);
            end
            for j=1:vectq(k)+1
                
                for t=p+q+1:n
                    if j==1
                        ddbeta(i,j,k)=ddbeta(i,j,k)+((tau(t,k)/(2*h(t,k)^2))*e(t-i+1,k)^2)/(n-p-q);
                    else
                        ddbeta(i,j,k)=ddbeta(i,j,k)+((tau(t,k)/(2*h(t,k)^2))*e(t-j+1,k)^2*e(t-i+1,k)^2)/(n-p-q);
                    end
                end
            end
        end
    end
    ddbeta(:,:,k)=inv(ddbeta(:,:,k));
end
for k=1:K
    for i=1:vectq(k)+1
        for j=1:vectq(k)+1
          BETA(i,k)=betaa(i,k)+ddbeta(i,j,k)*dbeta(i,k);  
        end
    end
end
% bet0=BETA(1,:);

for k=1:K
     if BETA(1,k)>0
        bet0(k)=BETA(1,k);
    else
        bet0(k)=-BETA(1,k);
    end
%     for i=1:vectq(k)
%         bet(i,k)=BETA(i+1,k);
%     end
     for i=1:vectq(k)
        if BETA(i+1,k)>0
            bet(i,k)=BETA(i+1,k);
        else
            bet(i,k)=-BETA(i+1,k);
        end
    end
end
end
