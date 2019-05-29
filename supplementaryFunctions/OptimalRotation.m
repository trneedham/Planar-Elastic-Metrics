function [p2new,angnew]=OptimalRotation(p1,p2,b)
    
    A = p1*p2';
    [U,~,V] = svd(A);
    if det(A)> 0
        Ot = U*V';
    else
        Ot = U*([V(:,1) -V(:,2)])';
    end
    angnew=atan2(Ot(1,2),Ot(1,1))*2*b;
    p2new = Ot*p2;
