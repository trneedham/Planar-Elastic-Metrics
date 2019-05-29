function p = FbTransform(c,b)

% Transforms a plane curve c into another plane curve p according to the
% b-transform. 

[~,n]=size(c);

v=curveDeriv(c);
[theta,rho,~]=curve2Polar(v);
psi=theta/(2*b);
nu=2*b*sqrt(rho);

[x(:),y(:)]=pol2cart(psi,nu);

p=zeros(2,n);

for j=1:n
    p(1,j)=x(j);
    p(2,j)=y(j);
end
