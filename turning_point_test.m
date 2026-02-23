function [H,PVAL]=turningpointtest(x,alpha)

n=length(x);
X(1:n-2,3)=x(1:n-2);
X(1:n-2,2)=x(2:n-1);
X(1:n-2,1)=x(3:n);

STATISTIC=sum((((X(:,2) > X(:,1)))&((X(:,2) > X(:,3))))|((((X(:,2) < X(:,1)) & ((X(:,2) < X(:,3)))))));

mu = 2*(n-2)/3;
sigma2 = (16*n-29)/90;
PVAL=2*(1-normcdf(abs(STATISTIC-mu)/sqrt(sigma2),0,1));
pz=norminv(1-alpha,0,1); 
H=abs(STATISTIC-mu)/sqrt(sigma2)>pz; %% 

