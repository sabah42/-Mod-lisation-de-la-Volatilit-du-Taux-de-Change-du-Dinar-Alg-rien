function [lfct] = opt(vect)
load donné x p q K n
%alph=vect(1:K)';
%phi0=vect(K+1:2*K)';
phi0=vect(1:K);
%ii=2*K;
ii=K;
for k=1:K
    for i=1:p
        phi(i,k)=vect(ii+i);
    end
    ii=ii+p;
end
beta0=vect(ii+1:ii+K)';
for k=1:K
    for i=1:q
        beta(i,k)=vect(ii+i+K);
    end
    ii=ii+q;
end

test2=1;
k=1;
while ( k<=K  && beta0(k)> 0 )
    k=k+1;
end
if k<K
    test2=0;
end
test3=1;
k=1;
while (k<=K && test3==1)
    i=1;
    while (i<=q && beta(i,k)>= 0)
        i=i+1;
    end
    if i~=q
        test3=0;
        k=k+1;
    end
end
if k~=K
    test3=0;
end
Test=test2*test3
if Test==1
    lfct=fct_vrai(alph,phi0,phi,beta0,beta);
else
    lfct=10^10;
end
end
