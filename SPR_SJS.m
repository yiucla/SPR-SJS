function [S,beta_S] = SPR_SJS(X,Y,delta,beta_ini,k)

iter = 1;
maxiter = 100;
tol = 1;
beta_old = beta_ini;
[n,p]= size(X);
u = Variance2(X,Y,delta,beta_ini);
gamma = zeros(p,1);

while (iter<=maxiter) && (tol>1E-3)
   
    m = X*beta_old;
    h = n^(-1/2);

    for i =1:p
        BB =repmat(delta',n,1).*(repmat(X(:,i),1,n)-repmat(X(:,i)',n,1)).*(repmat(Y,1,n)>=repmat(Y',n,1)).*normpdf((repmat(m,1,n)-repmat(m',n,1))/h);
        sb(i)=sum(sum(BB-diag(diag(BB))))/(n*(n-1)*h);
    end
    
    B = repmat(delta',n,1).*(repmat(Y,1,n)>=repmat(Y',n,1)).*normcdf((repmat(m,1,n)-repmat(m',n,1))/h);    
    lold = sum(sum(B-diag(diag(B))))/(n*(n-1));        
    gamma = beta_old + 1/u*sb';
    tmp = sort(abs(gamma),'descend');
    beta_new = gamma.*(abs(gamma)>=tmp(k));
    
    m = X*beta_new;
    C = repmat(delta',n,1).*(repmat(Y,1,n)>=repmat(Y',n,1)).*normcdf((repmat(m,1,n)-repmat(m',n,1))/h);  
    lnew = sum(sum(C-diag(diag(C))))/(n*(n-1));
    
    if lnew < lold
        u = 2*u;
    end
    tol = norm(beta_new-beta_old);
    beta_old = beta_new;
    iter = iter + 1;
end
S = find(beta_old~=0);
beta_S = beta_old(S);

return