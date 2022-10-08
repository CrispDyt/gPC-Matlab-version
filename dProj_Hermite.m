function y = dProj_Hermite(X, order, M, q)
%
% This routine projects a random variable in numeric form (discrete) onto
% Hermite-chaos basis.
%
% y = dProj_Hermite(X, order, M, q)
%
% Input:  X - random vector (row or column) of length N>>1;
%         order - Highest order of Hermite-chaos;
%         M - number of subdomains in random space (M=1 means no decomposition),
%             subdomains are equally spaced;
%         q - number of quadratures used to evaluate integrals (q = max(q,10));
% Output: 
%         y - coefficients of Hermite expansion of size M x (order+1), 
%         where each row in the expansion in each subdomain.
%
% Note:   if mod(N,M) ~= 0, then (M-1) points will be left out.
%
% Created by Dongbin Xiu 04/28/2003.
%

N=length(X);
Ni=floor(N/M);

q=max(q,10);
y=zeros(M,order+1);

SX=sort(X);
% h-discretization of samples
for m=1:M
    y(m,:)=dProj_Hermite1(SX((m-1)*Ni+1:m*Ni), order, q);
end
