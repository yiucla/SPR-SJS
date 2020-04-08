function u = Variance2(X,Y,delta,beta)

    [n,p] = size(X);
    A = zeros(p,p);
    m = X*beta;
    h = n^(-1/2);
    for i = 1:n
        for j=1:n
            if j ~= i
            A = A + delta(j)*(X(i,:)-X(j,:))'*(X(i,:)-X(j,:))*(Y(i)>=Y(j))*(m(i)-m(j))*normpdf((m(i)-m(j))/h)/(n*(n-1)*h^3);
            end
        end
    end
    u = max(eig(A));
    
    
return