function [phi,r,angles]=curve2Polar(v)

% Input: (samples of) parametric curve v
% Output: polar coordinate representation: phi is the true angle function,
% r is the radius function, angles is the list of exterior angles between
% the PL curve.

[~,n]=size(v);

[phiInitial,r]=cart2pol(v(1,:),v(2,:));


for j=1:n-1
    angles(j)=real(acos(dot(v(:,j)/norm(v(:,j)),v(:,j+1)/norm(v(:,j+1)))));
end

for j=1:n-1
    crossProds(:,j)=cross([v(:,j);0],[v(:,j+1);0]);
end

for j=1:n-1
    orientations(j)=sign(crossProds(3,j));
end

phi(1)=phiInitial(1);

for j=2:n
    phi(j)=phi(j-1)+orientations(j-1)*angles(j-1);
end
