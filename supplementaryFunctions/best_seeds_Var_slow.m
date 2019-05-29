function [theta2new,p2new,gamInew]=best_seeds_Var_slow(b,p1,p2,theta2)

% Input: A pair of planar curves p1 and p2 and a number k.
% Output: seeds = A list of the k best shifts of p2 for minimizing geodesic 
%         distance. The geodesic distance in this step is optimized over 
%         rotations, but not reparameterizations.

[~,T] = size(p2);
scl=8;
minE=1000;

for ctr = 0:floor((T)/scl)
    p2n = ShiftF(p2,ctr);
    theta2n = ShiftF_th(theta2,ctr);
    
    [gamI] = DynamicProgrammingQ(p2n/sqrt(InnerProd_Q(p2n,p2n)),p1/sqrt(InnerProd_Q(p1,p1)),0,0);
    gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    p2n = Group_Action_by_Gamma_Coordv2(p2n,gamI);
    theta2n = interp1(linspace(0,1,T) , theta2n ,gamI,'linear');
    
    [p2n,Ot]=OptimalRotation(p1,p2n,b);
    theta2n=theta2n-Ot;

    sc=sqrt(trapz(linspace(0,1,T),sum(p1.^2,1)));

    p1nn=p1/sc;
    p2nn=p2n/sc;

    dist=2*b*acos(trapz(linspace(0,1,T),sum(p1nn.*p2nn,1)));

    Ec=dist/(2*b);
    if Ec < minE
        p2new = p2n;
        theta2new=theta2n;
        minE = Ec;
        gamInew = gamI;
    end

end