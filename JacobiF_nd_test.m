function e2 = JacobiF_nd_test(ndim, norder, alpha, beta)
%
% JacobiF_nd.m - test routine JacobiF_nd by evaluating 
%                orthogonality via cubature rule.
%

ndim=3; norder=4;

[z,w] = ZWsmolyak_load(ndim,norder);
w=w/(2^ndim);

[npt,ntmp] = size(z);

P = nchoosek(ndim+norder,ndim);
poly=zeros(npt,P);

for m=1:npt
    poly(m,:) = JacobiF_nd(z(m,:),ndim,norder,0,0);
end

e2=zeros(1,P);
for m=1:P
    e2(m) = sum(poly(:,m).^2.*w);
end