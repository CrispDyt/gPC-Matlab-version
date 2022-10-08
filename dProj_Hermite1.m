function y = dProj_Hermite1(X, n, q)
%
%  function y = dProj_Hermite1(X, n, q)
%  
%  This function projects a random sample X onto the Hermite-chaos basis
%
%  Input:  X - random sample vector of certain length (N>>1),
%          n - highest-order of Hermite-chaos basis,
%          q - number of quadratures used for the projection. q=max(q,10);
%
%  Output: y - row vector of length (n+1) containing the coeffcients of
%              Hermite-chaos expansion, y(1)=y_0, and so on.
%
%  Code created by Dongbin Xiu on April 1, 2003.
%  Slight modification is made on 10/23/2003 to save the overwriting
%  and restoring of data X by using a dummy data Z.
%
q=max(q,10);

N=length(X);
y=zeros(1,n+1);

y(1)=mean(X);
sigma=std(X);
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
iF=norminv(s,0,1);

for k=2:n+1
y(k)=sum(iG.*HermiteF(iF,k-1).*w)/(2*factorial(k-1));
end
y(2:n+1)=y(2:n+1)*sigma;
