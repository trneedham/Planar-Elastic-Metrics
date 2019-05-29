function fn = Group_Action_by_Gamma_Coordv2(f,gamma)

[n,T] = size(f);

gammaGrad = gradient(gamma,1/T);

% for j=1:n
%     fn(j,:) = sqrt(gammaGrad).*interp1(linspace(0,1,T) , f(j,:) ,gamma,'linear');
% end

for j=1:n
    fn(j,:) = sqrt(gammaGrad).*spline(linspace(0,1,T) , f(j,:) ,gamma);
end