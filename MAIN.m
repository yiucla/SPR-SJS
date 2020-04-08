clc, clear
load('matlab.mat')
[n,p]=size(Z);
Y = X;
X = Z;
X = (X - ones(n,1)*mean(X))./(ones(n,1)*std(X));
Y = log(Y);
Y = Y - mean(Y);
delta = de;
d = floor(1.5*log(n)*n^(1/3));
mr = 1 - sum(delta)/n;

beta_ini = boosting(X,Y,delta);
[S,beta_S] = SPR_SJS(X,Y,delta,beta_ini,k); 

