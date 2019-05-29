function [dist,rho_geod,theta1,theta2i,theta2n,transform_geod,c_geod,gamI]=geodesic_dist(b,c1,c2,num_samp,num_steps,figures)

n=num_samp; %Number of samples in curve.
lam=0; %Parameter for Dynamic programming --- ignore.

%Resample, rescale, translate curves.
c1=preprocess_curve(ReSampleCurve(c1,n));
c2=preprocess_curve(ReSampleCurve(c2,n));

%Align the curves.
A = c1*c2';  
[U,~,V] = svd(A);
    if det(A)> 0
        Ot = U*V';
    else
        Ot = U*([V(:,1) -V(:,2)])';
    end
c2 = Ot*c2;

%Apply Fb Transform.
v1=curveDeriv(c1); %Find derivative curves.
v2=curveDeriv(c2);

[theta1,rho1,~]=curve2Polar(v1); %Send derivatives to polar coordinates.
[theta2,rho2,~]=curve2Polar(v2);

% Correct one of the angle functions to account for sharp angles
theta2i=theta2;
theta2new=theta2;

for j=1:n
    if theta2(j)-theta1(j)>2*pi*0.9
        for k=j:n
            theta2new(k)=theta2(k)-2*pi;
        end
    elseif theta2(j)-theta1(j)<-2*pi*0.9
        for k=j:n
            theta2new(k)=theta2(k)+2*pi;
        end    
    end
end

theta2=theta2new;

psi1=theta1/(2*b); %Apply transformation to new polar coordinates.
nu1=2*b*sqrt(rho1);

psi2=theta2/(2*b);
nu2=2*b*sqrt(rho2);

[x1(:),y1(:)]=pol2cart(psi1,nu1); %Send the result to rectangular coords.
[x2(:),y2(:)]=pol2cart(psi2,nu2);

p1=zeros(2,n);
p2=zeros(2,n);

for j=1:n
    p1(1,j)=x1(j);
    p1(2,j)=y1(j);
    p2(1,j)=x2(j);
    p2(2,j)=y2(j);
end

%Find optimal reparameterization.
[gamI] = DynamicProgrammingQ(p2/sqrt(InnerProd_Q(p2,p2)),p1/sqrt(InnerProd_Q(p1,p1)),lam,0);
gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
p2 = Group_Action_by_Gamma_Coordv2(p2,gamI);
theta2 = interp1(linspace(0,1,n) , theta2 ,gamI,'linear');
    
  
% Rotate the curves again.
[p2,Ot]=OptimalRotation(p1,p2,b);
theta2n=theta2-Ot;

% Compute the geodesic in transform space.
[dist,transform_geod]=Geod_Planar(p1,p2,b,num_steps);

% Compute polar representation of transform geodesic.
for k=1:num_steps
    [psi_geod(:,k),nu_geod(:,k),~]=curve2Polar(transform_geod(:,:,k));
end

% Send to geodesic in curve space. If the curves are complicated this may
% look bad.
for k=1:num_steps
    theta_geod(:,k)=2*b*psi_geod(:,k);
end

for k=1:num_steps
    rho_geod(:,k)=(nu_geod(:,k)/(2*b)).^2;
end

v_geod=zeros(2,n,num_steps);

for j=1:n
    for k=1:num_steps
        [v_geod(1,j,k),v_geod(2,j,k)]=pol2cart(theta_geod(j,k),rho_geod(j,k));
    end
end

c_geod=zeros(2,n,num_steps);

for k=1:num_steps
    for d=1:2
        c_geod(d,:,k)=preprocess_curve(cumsum(v_geod(d,:,k)));
    end
end        

% Displays reparameterization and geodesic, if figures option is selected.
if figures
    figure(1)
    clf
    plot(gamI,'LineWidth',3)
    title('Optimal Reparameterization')
    axis square
    for k=1:num_steps
        figure(2)
        clf
        plot(c_geod(1,:,k),c_geod(2,:,k))
        title('Geodesic (click through steps)')
        pause
    end
end