function c = FbTransformInverse(p,b)

% Takes a curve p in b-transform space to a plane curve c. For small values
% of b (<0.5 theoretically, but <0.4 in practice), this is not really an
% inverse to b because the angle function can be distorted.

[~,n] = size(p);

[psi,nu,~]=curve2Polar(p);

for k=1:n
    theta(k)=2*b*psi(k);
    rho(k)=(nu(k)/(2*b))^2;
end

for k=1:n
    [v(1,k),v(2,k)]=pol2cart(theta(k),rho(k));
end

for d=1:2
    c(d,:)=cumsum(v(d,:));
end
    