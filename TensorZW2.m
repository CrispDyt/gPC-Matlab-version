function [z,w] = TensorZW2(z1,w1, z2, w2)
%
% The routine takes TWO multi-d qudrature rules and makes their tensor
% product.
%
% z = TensorGridZW(z1, w1, z2, w2)
% 
% z1,w1: first qudrature rule. Each row is a point.
% z2,w2: second quadrature rule. Each row is a point.
%
% by Dongbin Xiu, 11/18/2013.
%

[nr1, nc1] = size(z1);  
[nr2, nc2] = size(z2); 


z = [z1(1,:) z2(1,:)];  w = w1(1)*w2(1);

for i=1:nr1
    for j=1:nr2
        if  ((i~=1) | (j~=1))
            z = [z; [z1(i,:) z2(j,:)]];
            w = [w; w1(i)*w2(j)];
        end
    end
end
