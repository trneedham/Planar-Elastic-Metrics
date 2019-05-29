%% 

% This uses the fast version of the program to create a geodesic. 
% As the b parameter is made small, the results will degrade.

clear all;

addpath('supplementaryFunctions')

load('comp_curves.mat')

b=0.4;

c1=C(:,:,1);
c2=C(:,:,3);

main_large_b(c1,c2,b)
%% 

% This uses the slower code which works for small values of b, but which 
% uses a gradient descent-based path straightening algorithm to produce the
% geodesic. 

clear all;

addpath('supplementaryFunctions')

load('comp_curves.mat')

b=0.2;

c1=C(:,:,1);
c2=C(:,:,3);

main_small_b(c1,c2,b)