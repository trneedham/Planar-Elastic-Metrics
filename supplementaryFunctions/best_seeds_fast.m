function seed_list = best_seeds_fast(c1,c2,b,num_seeds)

% Input: plane curves c1 and c2, metric parameter b, number of best seeds
% to take

% Output: list containing best seeds (of length num_seeds) without any
% reparameterization---just L2 dist in the transform space

[~,n] = size(c1);

% The first step perturbs the last point of the curves--if the first and
% last point are equal then this causes problems with the transforms.
if c1(:,n)==c1(:,1)
        c1(:,n)=7/8*c1(:,n)+1/8*c1(:,n-1);
end
    
if c2(:,n)==c2(:,1)
        c2(:,n)=7/8*c2(:,n)+1/8*c2(:,n-1);
end
    
p1 = FbTransform(c1,b);


dists = zeros(1,n);

for j=1:n
    c2Shifted = ShiftF(c2,j-1);
    p2 = FbTransform(c2Shifted,b);
    [p2,~]=OptimalRotation(p1,p2,b);
    [dist,~]=Geod_Planar(p1,p2,b,10);
    dists(j) = dist;
end

[~,Inds] = sort(dists);

for j=1:n
    Inds_fixed(j) = Inds(j)-1;
end

seed_list = Inds_fixed(1:num_seeds);
    