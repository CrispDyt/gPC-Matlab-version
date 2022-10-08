function pn = index_p1_hc( p, k, order )
% 
% This routine takes a row vector p and generates p1 whose entries
% are all the combination of those entries starting from k of p plus one.
%

d=length(p);
nrow = d-k+1;
p1 = zeros(nrow,d);
ptmp = p;

row=0;
for i=k:d
  ptmp = p;
  ptmp(i) = p(i) + 1;
  if norm_nix(ptmp) <= order
    row = row + 1;
    p1(row,:) = ptmp;
  end
end

if row == 0
    pn = [];
else
    pn = p1(1:row,:);
end
%p1 = sortrows(p1);
