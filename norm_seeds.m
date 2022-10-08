function  seeds = norm_seeds(X)
%
%  function seeds = norm_seeds(X)
%  
%  This function returns N(0,1) sequence which has the same seeds
%  as the input random sequence X.
%
%  Input:  X - random sample vector of certain length (N>>1),
%
%  Output: seeds - N(0,1) seeds with the same size of X.
%
%  Written by Dongbin Xiu on 10/23/2003.
%

N=length(X);   seeds = zeros(size(X));

y(1)=mean(X);
sigma=std(X);
if sigma < eps
   seeds = X-y(1);
   return;
end

%Z=(X-y(1))/sigma;       % normalize data
%Z=X;

% construct CDF of Z
[x,ix]=sort(X);
G=((1:N)-0.5)/N;

ss=norminv(G,0,1);

for i=1:N
  seeds(ix(i)) = ss(i);
end
