function [Phi0,Phi,Beta0,Beta,Alpha]=em2_ARCH(phi0,phi,alph,beta0,beta,x,vectp,vectq,I)

for i=1:I
    if i==1    
        [alp(i,:),ph0(i,:),Ph(:,:,i),bet0(i,:),bet(:,:,i)] =EM_step(phi0,phi,alph,beta0,beta,x,vectp,vectq);
        
    else
        [alp(i,:),ph0(i,:),Ph(:,:,i),bet0(i,:),bet(:,:,i)] =EM_step(ph0(i-1,:),Ph(:,:,i-1),alp(i-1,:),bet0(i-1,:),bet(:,:,i-1),x,vectp,vectq);
      
    end 
end
Phi0=ph0(end,:);
Phi=Ph(:,:,end);
Beta=bet(:,:,end);
Beta0=bet0(end,:);
Alpha=alp(end,:);
end