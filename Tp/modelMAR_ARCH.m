function [Alpha,Phi0,Phi,Beta0,Beta,BIC,H,pvalue,Qstat,CriticalValue,M]=modelMAR_ARCH(alph,phi0,phi,beta0,beta,x,vectp,vectq,I)

%[test1]=Stat_MAR(alph,phi);
%if test1==1
[Phi0,Phi,Beta0,Beta,Alpha]=em2_ARCH(phi0,phi,alph,beta0,beta,x,vectp,vectq,I);

%test2=Stat_MAR(Alpha,Phi);
% if test2==1
K=size(Alpha,2);
p=size(Phi,1);
q=size(Beta,1);
[phi0th,phith,beta0th,betath,alphth]=fish_ARCH(Alpha,Phi0,Phi,Beta0,Beta,x,vectp,vectq);

for k=1:K-1
    if abs(Alpha(k))*inv(alphth(k))<1.96
        display 'alpha est significativement égale à zéro'
        Alpha=0;
        
    else display 'alpha est significativement différent à zéro'
    end
end
for k=1:K
    if abs(Phi0(k)*inv(phi0th(k)))<1.96
        display 'phi0 est significativement égale à zéro'
        Phi0(k)=0;
    else  display 'phi0 est significativement différent à zéro'
    end
end
for k=1:K
    for j=1:vectp(k)
        if abs(Phi(j,k)*inv(phith(j,k)))<1.96
            display 'phi est significativement égale à zéro'
            Phi(j,k)=0;
        else  display 'phi est significativement différent à zéro'
        end
    end
end
for k=1:K
    if abs(Beta0(k)*inv(beta0th(k)))<1.96
        display 'beta0 est significativement égale à zéro'
        Phi0(k)=0;
    else  display 'beta0 est significativement différent à zéro'
    end
end
for k=1:K
    for j=1:vectq(k)
        if abs(Beta(j,k)*inv(betath(j,k)))<1.96
            display 'beta est significativement égale à zéro'
            Phi(j,k)=0;
        else  display 'beta est significativement différent à zéro'
        end
    end
end
[x_ajust,eps_ajust]=ajust_MAR_ARCH(x,Alpha,Phi0,Phi,Beta0,Beta,vectp,vectq);
[para,BIC]=critere_ARCH(Phi0,Phi,Beta0,Beta,Alpha,x,vectp,vectq);
[var]=var_MAR_ARCH(Alpha,Phi0,Phi,Beta0,Beta,x,vectp,vectq);

%************************ test sur les résidus ****************************

n=length(x);

figure(1)
fontSize=10;
subplot(3,1,1);
plot(x(p+q+1:n))
hold on
plot(x_ajust(p+q+1:n),'r-')
title(' Série ajusté ','FontSize', fontSize);
legend('Série','Série ajusté')
hold off
subplot(3,1,2);
plot(var(p+1:n))
title('Volatilité','FontSize', fontSize);
subplot(3,1,3);
plot(eps_ajust(p+q+1:n))
title('Série résiduelle','FontSize', fontSize);
figure(2)
subplot(2,1,1);
autocorr(eps_ajust(p+q+1:n))
subplot(2,1,2);
parcorr(eps_ajust(p+q+1:n))
figure(3)
 subplot(2,1,1);
 autocorr(eps_ajust(p+q+1:n).^2)
 subplot(2,1,2);
  parcorr(eps_ajust(p+q+1:n).^2)

%********************** nullité de la moyenne *****************************

M=(mean(eps_ajust(p+q+1:n))*sqrt(n-p-q));%  vérifier les paramètres non nul
if abs (M)<1.96
    display 'Le test de la nullité de la moyenne est vérifié'
else
    display 'Le test de la nullité de la moyenne n''est pas vérifié'
end

%************************** Test de ljung box *****************************

[H,pvalue,Qstat,CriticalValue]=lbqtest(eps_ajust(p+q+1:n),200,.05,200);
if H==0
    display 'l''hyothèse de non corrélation est acceptée'
else
    display 'l''hyothèse de non corrélation est rejetée'
end


end
%end
%end

