function y = dProj_Legendre1(X, n, q)
%
%  function y = dProj_Legendre1(X, n, q)
%  
%  This function projects a random sample X onto the Legendre-chaos basis
%
%  Input:  X - random sample vector of certain length (N>>1),
%          n - highest-order of Legendre-chaos basis,
%          q - number of quadratures used for the projection. q=max(q,10);
%
%  Output: y - row vector of length (n+1) containing the coeffcients of
%              Legendre-chaos expansion, y(1)=y_0, and so on.
%
%  Code created by Dongbin Xiu on September 15, 2004.
%
q=max(q,10);

N=length(X);
y=zeros(1,n+1);

y(1)=mean(X);
sigma=(max(X-y(1))-min(X-y(1)))/2;
if sigma < eps
   y(2)=sigma;
   return;
end

sigma=sigma*1.1;

Z=(X-y(1))/sigma;       % normalize data

% construct CDF of Z
x=sort(Z);
G=((1:N)-0.5)/N;

[z,w]=JacobiZW(q,0,0);
s=(1+z)/2;

iG=interp1(G,x,s,'cubic');
%iF=z;
J2 = jacobi_e2_1d(n, 0,0);

for k=2:n+1
y(k)=sum(iG.*JacobiF(z,k-1,0,0).*w)/(2*J2(k));
end
y(2:n+1)=y(2:n+1)*sigma;
