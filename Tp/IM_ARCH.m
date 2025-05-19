function [alph_vrai_v,alph_moy,alph_e_t_empirique,alph_e_t_theorique,phi_vrai_v,phi_moy,phi_e_t_empirique,phi_e_t_theorique,phi0_vrai_v,phi0_moy,phi0_e_t_empirique,phi0_e_t_theorique,beta0_vrai_v,beta0_moy,beta0_e_t_empirique,beta0_e_t_theorique,beta_vrai_v,beta_moy,beta_e_t_empirique,beta_e_t_theorique]=IM_ARCH(alph,phi0,phi,beta0,beta,vectp,vectq,nbr,n,I)
[p,K]=size(phi);
q=size(beta,1);
alph_vrai_v=alph;
phi_vrai_v=phi;
phi0_vrai_v=phi0;
beta_vrai_v=beta;
beta0_vrai_v=beta0;
for i=1:nbr
    i
    x(i,:) = MAR_ARCH( phi0,phi,alph,beta0,beta,n);
    [Phi0(i,:),Phi(:,:,i),Beta0(i,:),Beta(:,:,i),Alpha(i,:)]=em2_ARCH(phi0,phi,alph,beta0,beta,x(i,:),vectp,vectq,I);
    [phi0th(i,:),phith(:,:,i),beta0th(i,:),betath(:,:,i),alphth(i,:),dd]=fish_ARCH(Alpha(i,:),Phi0(i,:),Phi(:,:,i),Beta0(i,:),Beta(:,:,i),x(i,:),vectp,vectq); 
  bb=true;

end
for k=1:K
    alph_moy(k)=mean(Alpha(:,k));
    alph_e_t_empirique(k)=std(Alpha(:,k));
    phi0_moy(k)=mean(Phi0(:,k));
    phi0_e_t_empirique(k)=std(Phi0(:,k));
    beta0_moy(k)=mean(Beta0(:,k));
    beta0_e_t_empirique(k) =std(Beta0(:,k));
    for ii=1:p
        phi_moy(ii,k)=mean(Phi(ii,k,:));
        phi_e_t_empirique(ii,k)=std(Phi(ii,k,:));
    end
    for ii=1:q
        beta_moy(ii,k)=mean(Beta(ii,k,:));
        beta_e_t_empirique(ii,k)=std(Beta(ii,k,:));
    end
end
alph_e_t_theorique=alphth(end,:);
phi0_e_t_theorique=phi0th(end,:);
phi_e_t_theorique=phith(:,:,end);
beta0_e_t_theorique=beta0th(end,:);
beta_e_t_theorique=betath(:,:,end);

end
