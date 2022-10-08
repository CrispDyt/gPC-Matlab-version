function c = HermiteRing_1d(s, l)

c=zeros(s+l+1,1);

for r = abs(s-l):(s+l)
g=(l+s+r)/2;
if ( abs(g-round(g)) < 1e-3 )
c(r+1) = gamma(s+1)*gamma(l+1)/(gamma(g-l+1)*gamma(g-s+1)*gamma(g-r+1));
end
end


