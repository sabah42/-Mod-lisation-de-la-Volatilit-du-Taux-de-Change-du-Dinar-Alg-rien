function [Alpha,Phi0,Phi,Beta0,Beta,BIC,H,pvalue,Qstat,CriticalValue,M]=modelMARARCH(x,alpha,phi0,phi,beta0,beta,I)
%[test1,test2]=Stationarity_MARARCH(alpha,phi,beta);
%if test1==1 & test2==1

    [Alpha,Phi0,Phi,Beta0,Beta]=EM_MARARCH(x,alpha,phi0,phi,beta0,beta,itmax);
    %[test3,test4]=Stationarity_MARARCH(x,Alpha,Phi0,Phi,Beta0,Beta);
    if test3==1 & test4==1
        K=size(Alpha,2);
        p=size(Phi,1);
        q=size(Beta,1);
        for k=1:K-1
            if ebs(Alpha(k))*inv(sqrt(Diag(k)))<1.96
                display 'alpha est significativement égale à zéro'
                Alpha=0;
                
            else display 'alpha est significativement différent à zéro'
            end
        end
        for k=1:K
            if abs(Phi0(k)*inv(sqrt(Diag(K+(p+q+2)*(K-1)))))<1.96
                display 'phi0 est significativement égale à zéro'
                Phi0(k)=0;
            else  display 'phi0 est significativement différent à zéro'
            end
        end
        for k=1:K
            for j=1:p
                if abs(Phi(j,k)*inv(sqrt(Diag(K+j+(p+q+2)*(K-1)))))<1.96
                    display 'phi est significativement égale à zéro'
                    Phi(j,k)=0;
                else  display 'phi est significativement différent à zéro'
                end
            end
        end
        for k=1:K
            if abs(Beta0(k)*inv(sqrt(Diag(K+p+1+(p+q+2)*(K-1)))))<1.96
                display 'Beta0 est significativement égale à zéro'
                Beta0(k)=0;
            else  display 'Beta0 est significativement différent à zéro'
            end
        end
        for k=1:K
            for j=1:q
                if abs(Beta(j,k)*inv(sqrt(Diag(K+j+p+1+(p+q+2)*(K-1)))))<1.96
                    display 'Beta est significativement égale à zéro'
                    Beta(j,k)=0;
                else  display 'Beta est significativement différent à zéro'
                end
            end
        end
        [test1,test2]=Stationarity_MARARCH(Alpha,Phi,Beta);
        [x_ajust,eps_ajust]=ajust_MARARCH(x,Alpha,Phi0,Phi,Beta0,Beta);
        [BIC,para]=bic_MARARCH(x,Alpha,Phi0,Phi,Beta0,Beta);
        %%%%%%%%%%%test sur les résidus%%%%%%%%%%%%%%%%
        n=length(x);
        var2=var_marARCH(x,Alpha,Phi0,Phi,Beta0,Beta);
        figure(1)
        %parcorr(eps_ajust(p+q+1:n))
        plot(x(p+q+1:n))
        hold on
        plot(x_ajust(p+q+1:n),'r-')
        figure(2)
        plot(eps_ajust(p+q+1:n))
        figure(3)
        autocorr(eps_ajust(p+q+1:n))
        figure(4)
        autocorr(eps_ajust(p+q+1:n))
        figure(5)
        autocorr(eps_ajust(p+q+1:n).^2)
        figure(6)
        autocorr(eps_ajust(p+q+1:n).^2)
        figure(7)
        plot(var2(p+q+1:n))
        parcorr(eps_ajust(p+q+1:n),100)
        %%%%%%%%%%%%%%%%%%%%%%%%%%nullité de la moyenne%%%%%%%%%%%%%%%%%%%%%%%
        M=((mean(eps_ajust(p+q+1:n))*inv(std(eps_ajust(p+q+1:n)))*sqrt(n-para-p-q-1)));
        if abs ((mean(eps_ajust(p+q+1:n))*inv(std(eps_ajust(p+q+1:n)))*sqrt(n-para-p-q-1)))<1.96
            display 'Le test de la nullité de la moyenne est vérifié'
        else
            display 'Le test de la nullité de la moyenne n''est pas vérifié' 
        end
        %%%%%%%%%%Test de ljung box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [H,pValue,Qstat,CriticalValue]=lbqtest(eps_ajust(p+q+1:n),50,.05,50-para);
        if H==0
            display 'l''hyothèse de non corrélation est acceptée'
        else 
            display 'l''hyothèse de non corrélation est rejetée'
        end
    end
end
end
        
        
        
        
        
        
        
