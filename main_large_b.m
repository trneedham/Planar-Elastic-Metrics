function main_large_b(c1,c2,b)

num_samp=200; %Number of sample points in the curves
num_steps=7; %Number of steps to display in the geodesic
figures=1; %Show figures

[dist,~,~,~,~,~,~,~]=geodesic_dist(b,c1,c2,num_samp,num_steps,figures);

sprintf('The distance between the two curves is %0.3f',dist)