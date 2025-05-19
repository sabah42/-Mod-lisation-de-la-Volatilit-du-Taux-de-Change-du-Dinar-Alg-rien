function []=plo(x)


load('gbparch.mat')

fplot(@(y)fct_predic_arch( Phi0,Phi,Beta0,Beta,Alpha,2862,x,y),[-0.0001 0.0001])
load('gbpmar.mat')
phi0=Phi0;
phi=Phi;
sigma=Sigma;
alph=Alpha;
hold on
fplot(@(y)fct_predic( phi0,phi,sigma,alph,2862,x,y),[-0.0001 0.0001])
end
