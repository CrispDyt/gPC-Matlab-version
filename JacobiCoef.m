function coef = JacobiCoef( n, alpha, beta )
%
% JacobiCoef.m - returns the coefficients of Jacobi polynomials in monic form
%                P^n(x) = x^0 + a_1 x + ... + a_n x^n.
%                This is the format in matlab for polynomials.
%
% Syntax :   coef = JacobiCoef(n, alpha, beta)
%
% Input :    n - order of the Jacobi polynomial;
%            alpha, beta - input parameters defined by Jacobi polynomials.
% Output:    coef - (1 x (n+1)) vector containing a_i, the coefficients of 
%            x^(n-i), where i=0,...n. 
%
% by Dongbin Xiu 12/09/12.
%

coef = zeros(1,n+1);

for k=0:n
cc=zeros(1,n+1);
for m=0:k
cc(m+1)=(-1)^m/(factorial(m)*factorial(k-m));
%cc(m+1)=(-1)^m*nchoosek(k,m);
end
ak=pochhammer(n+alpha+beta+1,k)*pochhammer(alpha+k+1,n-k)*pochhammer(-n,k)/(2^k);
coef=coef+ak*cc;
end

coef=fliplr(coef/factorial(n));
