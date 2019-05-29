function [dist,geod]=Geod_Planar(p1,p2,b,k)

% Input: two planar curves in b-transfrom space and number of steps to
% record in the geodesic joining them.

% Output: geodesic distance in the Hilbert sphere and the geodesic between
% them.

[~,N]=size(p1);

sc=sqrt(trapz(linspace(0,1,N),sum(p1.^2,1)));

p1n=p1/sc;
p2n=p2/sc;

dist=2*b*acos(trapz(linspace(0,1,N),sum(p1n.*p2n,1)));

theta=dist/(2*b);

% Define the parameter for the geodesic path.

u=linspace(0,1,k);

for j=1:k
    geod(:,:,j)=(sc/sin(theta))*(sin((1-u(j))*theta)*p1n+sin(u(j)*theta)*p2n);
end