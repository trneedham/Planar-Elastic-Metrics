function [d,dists] = GeodesicElasticClosedVarParams(p1,p2,b)

% input p1 and p2 as 2xn or 3xn matrices (planar vs. 3D closed curves)
% to turn off figures set figs=0
% the output is the distance d between p1 and p2 and the geodesic path Geod
% stp refers to the number of shapes displayed along the geodesic


N = 100; % Resample the curves to have the same number of points.

p1 = preprocess_curve(ReSampleCurve(p1,N));
p2 = preprocess_curve(ReSampleCurve(p2,N));

q1 = FbTransform(p1,b);
q2 = FbTransform(p2,b);

[~,~,dists] = Find_Rotation_and_Seed_unique_v2(p1,p2,b);
%[~,~,dists,~] = Find_Rotation_and_Seed_unique(q1,q2,p2,b,1);

d = min(dists);



