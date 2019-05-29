function [q2best,Rbest,dists,gamIbest] = Find_Rotation_and_Seed_unique_v2(p1,p2,b)

[~,T] = size(p1);

if p1(:,T) == p1(:,1)
    p1(:,T) = 7/8*p1(:,T) + 1/8*p1(:,T-1);
end

if p2(:,T) == p2(:,1)
    p2(:,T) = 7/8*p2(:,T) + 1/8*p2(:,T-1);
end

q1 = FbTransform(p1,b);
q1 = q1/norm_Q(q1);

v2 = curveDeriv(p2);
[theta,rho,~]=curve2Polar(v2);


scl = 1;
minE = 1000;

for ctr = 0:floor((T)/scl)
    p2n = ShiftF(p2,scl*ctr);
    thetan = ShiftF(theta,scl*ctr);
    rhon = ShiftF(rho,scl*ctr);
    psin=thetan/(2*b);
    nun=2*b*sqrt(rhon);
    [xn(:),yn(:)]=pol2cart(psin,nun);
    q2n = [xn;yn];
%    q2n = FbTransform(p2n,b);
%    q2n = q2n/norm_Q(q2n);
%     q2 = FbTransform(p2,b);
%     q2 = q2/norm_Q(q2);
%    q2n = ShiftF(q2,scl*ctr);
    [q2n,R] = Find_Best_Rotation(q1,q2n);
    
        
        if norm(q1-q2n,'fro') > 0.0001
            gam = DynamicProgrammingQ(q1,q2n,0,0);
            gamI = invertGamma(gam);
            gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
            p2new = Group_Action_by_Gamma_Coord(p2n,gamI);
            v2new = curveDeriv(p2new);
            [~,rhonew,~]=curve2Polar(v2new);
            thetanew = Group_Action_by_Gamma_Coord(thetan,gamI);
            %rhonew = ShiftF(rhonew,scl*ctr);
            psinew=thetanew/(2*b);
            nunew=2*b*sqrt(rhonew);
            [xnew(:),ynew(:)]=pol2cart(psinew,nunew);
            q2new = [xnew;ynew];
            %q2new = R*FbTransform(p2new,b);
        else
            q2new = q2n;
        end

    
    Ec = acos(InnerProd_Q(q1/norm_Q(q1),q2new/norm_Q(q2new)));
    
    dists(ctr+1) = Ec;
    
    if Ec < minE
        Rbest=R;
        q2best  = q2new;
        minE = Ec;
        gamIbest = gamI;
    end
end