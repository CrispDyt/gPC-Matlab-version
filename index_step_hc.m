function pn = index_step_hc( p, order )
% 
% This routine takes a matrix p and generates matrix p1 whose entries
% are all the combination of those of rows of p plus 1 to n.
%

pn = p;
p1 = index_step1_hc(p, order);

while isempty(p1) ~= 1
   pn = cat(1, pn, p1);
   p1 = index_step1_hc(p1, order);
end