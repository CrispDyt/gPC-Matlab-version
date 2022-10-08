function f = pochhammer(a, n)
%
% pochhammer.m - returns value of the pochhammer symbol (a)_n, where n>=0.
%
% by Dongbin Xiu 12/09/02.
%

if n == 0
f = 1;
else
f=prod(a:a+n-1,2);
end
