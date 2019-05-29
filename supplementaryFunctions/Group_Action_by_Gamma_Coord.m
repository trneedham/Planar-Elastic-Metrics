function fn = Group_Action_by_Gamma_Coord(f,gamma)

[n,T] = size(f);

for j=1:n
    fn(j,:) = spline(linspace(0,1,T) , f(j,:) ,gamma);
end

return;