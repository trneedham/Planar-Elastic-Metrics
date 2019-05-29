function main_small_b(c1,c2,b)

num_samp=200; %Number of sample points in the curves
num_steps=7; %Number of steps to display in the geodesic
figures=0; %Leave as 0
maxiter=300; %Number of iterations in the path straightening algorithm. 

B=form_basis_path_straight(num_steps,num_samp,150,num_steps,num_steps);

[dist,rho_geod,theta1,~,theta2n,~,~,gamI]=geodesic_dist(b,c1,c2,num_samp,num_steps,figures);

sprintf('The distance between the two curves is %0.3f',dist)

figure(1)
clf
plot(gamI,'LineWidth',3)
title('Optimal Reparameterization')
axis square
    
[~,~,~,~] = path_straight(rho_geod',theta1,theta2n,B,maxiter);