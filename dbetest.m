function [dbeta,ddbeta] = dbetest( ph0,Ph,beta,beta0,tau,x,vectp,vectq)
[n,K]=size(tau);
p=size(Ph,1);
q=size(beta,1);
dbeta=zeros(q+1,K);
ddbeta=zeros(q+1,q+1,K);
for k=1:K
    for t=p+1:n
        e(t,k)=x(t)-ph0(k);
        for i=1:vectp(k)
            e(t,k)=e(t,k)-(Ph(i,k)*x(t-i));
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
    for t=p+q+1:n
        dhdbeta(1,t,k)=1;
        for i=1:vectq(k)
            dhdbeta(i+1,t,k)=e(t-i,k)^2;
        end
    end
end

for k=1:K
    for t=p+q+1:n
        dbeta(1,k)=dbeta(1,k)+(tau(t,k)/(2*h(t,k))*dhdbeta(1,t,k)*((e(t,k)^2/h(t,k))-1))/(n-p-q);
        ddbeta(1,1,k)=ddbeta(1,1,k)-(tau(t,k)/(2*h(t,k)^2)*(dhdbeta(1,t,k))^2)/(n-p-q);
        for i=1:vectq(k)
            dbeta(i+1,k)=dbeta(i+1,k)+(tau(t,k)/(2*h(t,k))*dhdbeta(i+1,t,k)*((e(t,k)^2/h(t,k))-1))/(n-p-q);
            ddbeta(i+1,1,k)=ddbeta(i+1,1,k)-(tau(t,k)/(2*h(t,k)^2)*dhdbeta(1,t,k)*dhdbeta(i+1,t,k))/(n-p-q);
            for j=1:vectq(k)
                ddbeta(i+1,j+1,k)=ddbeta(i+1,j+1,k)-(tau(t,k)/(2*h(t,k)^2)*dhdbeta(j+1,t,k)*dhdbeta(i+1,t,k))/(n-p-q);
            end
        end
        for j=1:vectq(k)
            ddbeta(1,j+1,k)=ddbeta(1,j+1,k)-(tau(t,k)/(2*h(t,k)^2)*dhdbeta(j+1,t,k)*dhdbeta(1,t,k))/(n-p-q);
        end   
    end
    
end
end
