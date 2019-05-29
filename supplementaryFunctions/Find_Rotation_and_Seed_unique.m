function [q2best,Rbest,dists,gamIbest] = Find_Rotation_and_Seed_unique(q1,q2,p2,b,reparamFlag)

[~,T] = size(q1);

scl = 1;
minE = 1000;
for ctr = 0:floor((T)/scl)
    p2n = ShiftF(p2,scl*ctr);
    q2n = FbTransform(p2n,b);
    %q2n = ShiftF(q2,scl*ctr);
    [q2n,R] = Find_Best_Rotation(q1,q2n);
    
    if(reparamFlag)
        
        if norm(q1-q2n,'fro') > 0.0001
            gam = DynamicProgrammingQ(q1,q2n,0,0);
            gamI = invertGamma(gam);
            gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
            p2new = Group_Action_by_Gamma_Coord(p2n,gamI);
            q2new = R*FbTransform(p2new,b);
        else
            q2new = q2n;
        end
        
    else
        q2new  = q2n;
    end
    
    Ec = acos(InnerProd_Q(q1,q2new));
    
    dists(ctr+1) = Ec;
    
    if Ec < minE
        Rbest=R;
        q2best  = q2new;
        minE = Ec;
        gamIbest = gamI;
    end
end