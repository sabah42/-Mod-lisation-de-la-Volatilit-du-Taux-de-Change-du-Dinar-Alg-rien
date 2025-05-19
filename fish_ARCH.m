function[phi0th,phith,beta0th,betath,alphth]=fish_ARCH(Alpha,Phi0,Phi,Beta0,Beta,x,vectp,vectq)
n=length(x);
K=size(Alpha,2);
p=size(Phi,1);
q=size(Beta,1);
IC0=zeros(K-1,K-1);
IM0=zeros(K-1,K-1);

%************************ calcul de l'erreur ******************************
for t=p+1:n
    for k=1:K
        e(t,k)=x(t)-Phi0(k);
        for i=1:vectp(k)
            e(t,k)=e(t,k)-Phi(i,k)*x(t-i);
        end
    end
end
%*********************** Espérance conditionnelle *************************
for t=p+q+1:n
    for k=1:K
        h(t,k)=Beta0(k);
        for j=1:vectq(k)
            h(t,k)=h(t,k)+Beta(j,k)*((e(t-j,k))^2);
        end
        S(k)=Alpha(k)*inv(sqrt(h(t,k)))*normpdf(e(t,k)*inv(sqrt(h(t,k))));
    end
    tau(t,:)=S/sum(S);
end
% ************************ les dérivées utilisés **************************
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
            E1(j)=Beta(j,k)*e(t-j,k);
        end
        dhdphi(1,t,k)=-2*sum(E1);
        for i=1:vectp(k)
            for j=1:vectq(k)
                E2(j)=Beta(j,k)*e(t-j,k)*dedphi(i+1,t-j,k);
            end
            dhdphi(i+1,t,k)=2*sum(E2);
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
%****************** la matrice IM0 && la matice IC0 ***********************
for l=1:K-1
    for k=1:K-1
        if (l==k)
            for t=p+q+1:n
                IC0(l,k)=IC0(l,k)+(tau(t,k)/Alpha(k)^2)+(tau(t,K)/Alpha(K)^2);
                IM0(l,k)=IM0(l,k)+(tau(t,k)*(1-tau(t,k))/Alpha(k)^2)+(tau(t,K)*((1-tau(t,K))/Alpha(K)^2))+(2*tau(t,k)*tau(t,K)/(Alpha(k)*Alpha(K)));
            end
        else
            
            for t=p+q+1:n
                IC0(l,k)=IC0(l,k)+(tau(t,K)/Alpha(K)^2);
                IM0(l,k)=IM0(l,k)+((tau(t,k)*tau(t,K)/(Alpha(k)*Alpha(K)))+(tau(t,l)*tau(t,K))/(Alpha(l)*Alpha(K))-((tau(t,l)*tau(t,k))/(Alpha(l)*Alpha(k)))+(tau(t,K)*(1-tau(t,K)))/(Alpha(K)^2));
            end
        end
    end
end
%************************** les blocs diagonaux ***************************

for k=1:K
    IC1=zeros(vectp(k)+1,vectp(k)+1,K);
    IC2=zeros(vectq(k)+1,vectq(k)+1,K);
    IMKTT=zeros(vectp(k)+1,vectp(k)+1,K);
    IMKBB=zeros(vectq(k)+1,vectq(k)+1,K);
    IMKBT=zeros(vectq(k)+1,vectp(k)+1,K);
    
end
for k=1:K
    for i=0:vectp(k)
        for j=0:vectp(k)
            for t=p+q+1:n
                IC1(i+1,j+1,k)=IC1(i+1,j+1,k)+(tau(t,k)*((dhdphi(i+1,t,k)*dhdphi(j+1,t,k)*inv(2*h(t,k)^2))+(dedphi(i+1,t,k)*dedphi(j+1,t,k)*inv(h(t,k)))));
                IMKTT(i+1,j+1,k)=IMKTT(i+1,j+1,k)+(tau(t,k)*(1-tau(t,k))*(dhdphi(i+1,t,k)*inv(2*h(t,k))*(e(t,k)^2*inv(h(t,k))-1)-(e(t,k)*inv(h(t,k)))*dedphi(i+1,t,k))*((dhdphi(j+1,t,k)*inv(2*h(t,k))*(e(t,k)^2*inv(h(t,k))-1)-(e(t,k)*inv(h(t,k))*dedphi(j+1,t,k)))));
                
            end
        end
    end
end
for k=1:K
    for i=0:vectq(k)
        for j=0:vectq(k)
            for t=p+q+1:n
                IC2(i+1,j+1,k)=IC2(i+1,j+1,k)+(tau(t,k)*inv(2*h(t,k)^2)*dhdbeta(i+1,t,k)*dhdbeta(j+1,t,k));
                IMKBB(i+1,j+1,k)=IMKBB(i+1,j+1,k)+(tau(t,k)*(1-tau(t,k))*inv(4*h(t,k)^2)*dhdbeta(i+1,t,k)*dhdbeta(j+1,t,k)*(e(t,k)^2*inv(h(t,k))-1)^2);
            end
        end
    end
end
for k=1:K
    for i=0:vectq(k)
        for j=0:vectp(k)
            for t=p+q+1:n
                IMKBT(i+1,j+1,k)=IMKBT(i+1,j+1,k)+(tau(t,k)*(1-tau(t,k))*inv(2*h(t,k))*dhdbeta(i+1,t,k)*(e(t,k)^2*inv(h(t,k))-1)*(dhdphi(j+1,t,k)*inv(2*h(t,k))*(e(t,k)^2*inv(h(t,k))-1)-(e(t,k)*inv(h(t,k))*dedphi(j+1,t,k))));
            end
            
        end
    end
end

for k=1:K
    IMKK=zeros(vectp(k)+vectp(k)+2,vectp(k)+vectp(k)+2,K) ;
end
%****************** remplir la matice IMKK ********************************
for k=1:K
    
    for i=1:vectq(k)+1
        for j=1:vectp(k)+1
            IMKK(i+vectp(k)+1,j,k)=IMKBT(i,j,k);
        end
    end
    IMKK(:,:,k)=IMKK(:,:,k)+IMKK(:,:,k)';
end
for k=1:K
    for i=1:vectp(k)+1
        for j=1:vectp(k)+1
            IMKK(i,j,k)=IMKTT(i,j,k);
        end
    end
    for i=1:vectq(k)+1
        for j=1:vectq(k)+1
            IMKK(i+vectp(k)+1,j+vectp(k)+1,k)=IMKBB(i,j,k);
        end
    end
    
end
%********************* remplir la matrice ICK (complet)********************
ICK=zeros(sum(vectp)+sum(vectq)+3*K-1);
for i=1:K-1
    for j=1:K-1
        ICK(i,j)=IC0(i,j);
    end
end
ii=K-1;
for k=1:K
    for i=1:vectp(k)+1
        for j=1:vectp(k)+1
            ICK(i+ii,j+ii)=IC1(i,j,k);
        end
    end
    ii=ii+vectp(k)+1;
    for i=1:vectq(k)+1
        for j=1:vectq(k)+1
            ICK(i+ii,j+ii)=IC2(i,j,k);
        end
    end
    ii=ii+vectq(k)+1;
end
IMM=zeros(sum(vectp)+sum(vectq)+3*K-1);

%************************* la matrice IMk0 ********************************
for k=1:K-1
    IM0T=zeros(vectp(k)+1,K-1,K-1);
end
for k=1:K-1
    for i=0:vectp(k)
        for t=p+q+1:n
            
            IM0T(i+1,k,k)=IM0T(i+1,k,k)+(((tau(t,k)*(1-tau(t,k))*inv(Alpha(k)))+tau(t,k)*tau(t,K)*inv(Alpha(K)))*(dhdphi(i+1,t,k)*inv(2*h(t,k))*((e(t,k)^2*inv(h(t,k)))-1)-(e(t,k)*inv(h(t,k))*dedphi(i+1,t,k))));
        end
        for L=1:K-1
            if L~=k
                for t=p+q+1:n
                    IM0T(i+1,L,k)=IM0T(i+1,L,k)+((((tau(t,k)*tau(t,K))*inv(Alpha(K)))-(tau(t,k)*tau(t,L)*inv(Alpha(L))))*(dhdphi(i+1,t,k)*inv(2*h(t,k))*((e(t,k)^2*inv(h(t,k))-1)-(e(t,k)*inv(h(t,k))*dedphi(i+1,t,k)))));
                end
                
            end
        end
    end
end
for k=1:K-1
    IM0B=zeros(vectq(k)+1,K-1,K-1);
end

for k=1:K-1
    for i=0:vectq(k)
        for t=p+q+1:n
            IM0B(i+1,k,k)=IM0B(i+1,k,k)+(((tau(t,k)*(1-tau(t,k))*inv(Alpha(k)))+(tau(t,k)*tau(t,K)*inv(Alpha(K))))*(inv(2*h(t,k))*dhdbeta(i+1,t,k)*((e(t,k)^2*inv(h(t,k)))-1)));
        end
        
        for L=1:K-1
            if L~=k
                
                for t=p+q+1:n
                    IM0B(i+1,L,k)=IM0B(i+1,L,k)+((((tau(t,k)*tau(t,K))*inv(Alpha(K)))-(tau(t,k)*tau(t,L)*inv(Alpha(L))))*(inv(2*h(t,k))*dhdbeta(i+1,t,k))*((e(t,k)^2*inv(h(t,k)))-1));
                end
                
            end
        end
    end
end

for k=1:K-1
    IMk0=zeros(vectp(k)+vectq(k)+2,K-1,K-1);
end
for k=1:K-1
    for j=1:K-1
        for i=1:vectp(k)+1
            IMk0(i,j,k)=IM0T(i,j,k);
        end
        for i=1:vectq(k)+1
            IMk0(i+vectp(k)+1,j,k)=IM0B(i,j,k);
        end
    end
end
%*************************** la matrice IMK0 ******************************
for k=1:K-1
    IM0KT=zeros(vectp(k)+1,K-1,K-1);
end
for k=1:K-1
    for i=0:vectp(k)
        
        for t=p+q+1:n
            
            IM0KT(i+1,k,k)=IM0KT(i+1,k,k)+(((-tau(t,K)*(1-tau(t,K))*inv(Alpha(K)))-(tau(t,k)*tau(t,K)*inv(Alpha(k))))*((dhdphi(i+1,t,K)*inv(2*h(t,K)))*((e(t,K)^2*inv(h(t,K))-1))-((e(t,K)*inv(h(t,K)))*dedphi(i+1,t,K))));
        end
        
    end
end
for k=1:K-1
    IM0KB=zeros(vectq(k)+1,K-1,K-1);
end
for k=1:K-1
    for i=0:vectq(k)
        for t=p+q+1:n
            IM0KB(i+1,k,k)=IM0KB(i+1,k,k)+(((-tau(t,K)*(1-tau(t,K))/Alpha(K))-(tau(t,k)*tau(t,K)/Alpha(k)))*(inv(2*h(t,K)))*dhdbeta(i+1,t,K)*((e(t,K)^2*inv(h(t,K))-1)));
        end
    end
end
%********************* affectation de la matrice IMK0 *********************
for k=1:K-1
    IMK0=zeros(vectp(k)+vectq(k)+2,K-1,K-1);
end
for k=1:K-1
    for j=1:K-1
        for i=1:vectp(k)+1
            IMK0(i,j,k)=IM0KT(i,j,k);
        end
        for i=1:vectq(k)+1
            IMK0(i+vectp(k)+1,j,k)=IM0KB(i,j,k);
        end
    end
end
jj=K-1;
for k=1:K-1
    for i=1:vectp(k)+vectq(k)+2
        IMM(i+jj,1:K-1)=IMk0(i,1:K-1,k);
    end
    jj=jj+vectp(k)+vectq(k)+2;
end
ii=sum(vectp)+sum(vectq)+3*(K-1)-vectp(K)-vectq(K);
for i=1:vectp(K)+vectq(k)+2
    IMM(i+ii,1:K-1)=IMK0(i,1:K-1);
end
for l=1:K
    for k=1:K-1
        IMKLTT=zeros(vectp(k)+1,vectp(l)+1,K,K);
    end
end
for l=1:K
    for k=2:K
        if k~=l
            for i=0:vectp(k)
                for j=0:vectp(l)
                    for t=p+q+1:n
                       IMKLTT(i+1,j+1,k,l)=IMKLTT(i+1,j+1,k,l)-((tau(t,k)*tau(t,l))*(((inv(2*h(t,k))*dhdphi(i+1,t,k))*((((e(t,k))^2)*inv(h(t,k)))-1))-(e(t,k)*inv(h(t,k)))*dedphi(i+1,t,k))*((dhdphi(j+1,t,l)*inv(2*h(t,l))*((e(t,l)^2*inv(h(t,l)))-1))-((e(t,l)*inv(h(t,l)))*dedphi(j+1,t,l))));   
                    end
                end
            end
        end
    end
end
for l=1:K
    for k=1:K-1
        IMKLBB=zeros(vectq(k)+1,vectq(l)+1,K,K);
    end
end
for l=1:K
    for k=2:K
        if k~=l
            for i=0:vectq(k)
                for j=0:vectq(l)
                    for t=p+q+1:n
                        IMKLBB(i+1,j+1,k,l)=IMKLBB(i+1,j+1,k,l)+(-(tau(t,k)*tau(t,l))*((dhdbeta(i+1,t,k)*inv(2*h(t,k)))*((e(t,k)^2*inv(h(t,k)))-1))*((inv(2*h(t,l)))*dhdbeta(j+1,t,l)*((e(t,l)^2*inv(h(t,l)))-1)));
                    end
                end
            end
        end
    end
end
for l=1:K
    for k=1:K-1
        IMKLBT=zeros(vectq(k)+1,vectp(l)+1,K,K);
    end
end
for l=1:K
    for k=2:K
        if k~=l
            for i=0:vectq(k)
                for j=0:vectp(l)
                    for t=p+q+1:n
                        IMKLBT(i+1,j+1,k,l)=IMKLBT(i+1,j+1,k,l)+(-tau(t,k)*tau(t,l)*(((inv(2*h(t,k)))*dhdbeta(i+1,t,k)*((e(t,k)^2*inv(h(t,k)))-1))*((inv(2*h(t,l))*dhdphi(j+1,t,l)*((e(t,l)^2*inv(h(t,l)))-1)-(e(t,l)*inv(h(t,l))*dedphi(j+1,t,l))))));
                    end
                end
            end
        end
    end
end



for l=1:K
    for k=1:K-1
        IMKLTB=zeros(vectp(k)+1,vectq(l)+1,K,K);
    end
end
for l=1:K
    for k=2:K
        if k~=l
            for i=0:vectp(k)
                for j=0:vectq(l)
                    for t=p+q+1:n
                        IMKLTB(i+1,j+1,k,l)=IMKLTB(i+1,j+1,k,l)+(-tau(t,k)*tau(t,l)*(((inv(2*h(t,k))*dhdphi(i+1,t,k))*((e(t,k)^2*inv(h(t,k)))-1)-(e(t,k)*inv(h(t,k)))*dedphi(i+1,t,k))*((inv(2*h(t,l)))*dhdbeta(j+1,t,l)*((e(t,l)^2*inv(h(t,l)))-1))));
                    end
                end
            end
        end
    end
end

%remplir la matrice IMKL
for l=1:K
    for k=2:K
        IMKL=zeros(vectq(k)+vectp(k)+2,vectp(l)+vectq(l)+2,K,K);
    end
end
for l=1:K
    for k=2:K
        for i=1:vectp(k)+1
            for j=1:vectp(k)+1
                IMKL(i,j,k,l)=IMKLTT(i,j,k,l);
            end
        end
        for i=1:vectq(k)+1
            for j=1:vectp(k)+1
                IMKL(i+vectp(k)+1,j,k,l)=IMKLBT(i,j,k,l);
            end
        end
        for i=1:vectq(k)+1
            for j=1:vectq(k)+1
                IMKL(i+vectp(k)+1,j+vectp(k)+1,k,l)=IMKLBB(i,j,k,l);
            end
        end
        for i=1:vectp(k)+1
            for j=1:vectq(k)+1
                IMKL(i,j+vectp(k)+1,k,l)=IMKLTB(i,j,k,l);
            end
        end
    end
    
end
jj=K-1;
for l=1:K-1
    ii=vectp(1)+vectq(1)+K+1;
    for k=2:K
        if k~=l
            for j=1:vectp(l)+vectq(l)+2
                for i=1:vectp(k)+vectq(k)+2
                    IMM(i+ii,j+jj)=IMKL(i,j,k,l);
                end
            end
        end
        ii=ii+vectp(k)+vectq(k)+2;
    end
    
    jj=jj+vectp(l)+vectq(l)+2;
end
IMM=IMM+IMM';
for i=1:K-1
    for j=1:K-1
        IMM(i,j)=IM0(i,j);
    end
end
ii=K-1;
for k=1:K
    for i=1:vectp(k)+vectq(k)+2
        for j=1:vectp(k)+vectq(k)+2
            IMM(i+ii,j+ii)=IMKK(i,j,k);
        end
    end
    ii=ii+vectp(k)+vectq(k)+2;
end
I=ICK-IMM;
I=inv(I);
dd=diag(I);
d=sqrt(dd);
%***************  détermination de l'ecart type théorique *****************
alphth=d(1:K-1);
alphth(K)=0;
for i=1:K-1
    for j=1:K-1
        alphth(K)=alphth(K)+I(i,j);
    end
end
alphth(K)=sqrt(alphth(K));
phi0th(1)=d(K);
beta0th(1)=d(K+vectp(1)+1);
j=K;
for k=2:K
    phith(1:vectp(k-1),k-1)=d(j+1:(j+vectp(k-1)));
    phi0th(k)=d(j+vectp(k-1)+vectq(k-1)+2);
    beta0th(k)=d(j+2*vectp(k-1)+vectq(k-1)+3);
    betath(1:vectq(k-1))=d(j+vectp(k-1)+2:(j+vectp(k-1)+vectq(k-1)+1));
    j=j+vectp(k-1)+vectq(k-1)+2;
end
phith(1:vectp(K),K)=d(j+1:(j+vectp(K)));
betath(1:vectq(K),K)=d(j+vectp(K)+2:j+vectp(K)+vectq(K)+1);
end
