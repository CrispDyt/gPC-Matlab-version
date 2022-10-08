function e = hermite_h3( im, jm, km )
%
% hermite_h3.m - Evaluate the inner product of Hermite-chaos basis
%                <H_im H_jm H_km>, where i,j,k are multi-dimensional indices
%                containing its order in each dimension.
%                The exact formula is employed.
%
% Syntax     e = hermite_h3( im, jm, km )
%
% Input:     im,jm,km = (1 x M) arrays with each entry denoting its order in
%                    that dimension.
%
% Output:    e = the result in integer form.
% 
% By Dongbin Xiu   2/23/2004.
%

M = length(im);

e = 1;
for m = 1:M
   i = im(m);  j = jm(m);  k = km(m);
   s = i + j + k;
   if mod(s,2) == 0
      s = s / 2;
      if ((s >= i) & (s>=j) & (s>=k))
	 em = factorial(i)*factorial(j)*factorial(k)/ ...
                 (factorial(s-i)*factorial(s-j)*factorial(s-k));
      else
         em = 0;
      end
   else
      em = 0; 
   end
   if em == 0
      e = 0;
      break;
   else
      e = e * em;
   end
end
