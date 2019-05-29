function [Thet,Thet_init,v_geod,v_geodinit,c_geod] = path_straight(R,th1,th2,B,maxiter)

S=size(R,1);

T=length(th1);

s=linspace(0,1,S);
t=linspace(0,1,T);

Nb=size(B,3);

for j=1:length(s)
    Thet(j,:)=(1-s(j))*th1+s(j)*(th2);
end

Thet_init=Thet;

[Thett,Thets]=gradient(Thet,1/(T-1),1/(S-1));

for j=1:Nb
    [~,Bs(:,:,j)]=gradient(B(:,:,j),1/(T-1),1/(S-1));
end

maxiter=300;

eps=0.05;

for k=1:maxiter

gradE=zeros(S,T);
    
for j=1:Nb
    dirderb=trapz(t,trapz(s,Thets.*Bs(:,:,j).*R))*B(:,:,j);
    gradE=gradE+dirderb;
    gradE(1,:)=zeros(1,T);
    gradE(end,:)=zeros(1,T);
end

Thet=Thet-eps*gradE;

[Thett,Thets]=gradient(Thet,1/(T-1),1/(S-1));

ngradE=sqrt(trapz(t,trapz(s,gradE.^2)));

nE=trapz(t,trapz(s,Thets.^2.*R));

[ngradE nE]

end

for j=1:S
    for k=1:T
        [v_geod(1,k,j),v_geod(2,k,j)]=pol2cart(Thet(j,k),R(j,k));
        [v_geodinit(1,k,j),v_geodinit(2,k,j)]=pol2cart(Thet_init(j,k),R(j,k));
    end
end

c_geod=zeros(2,T,S);

for k=1:S
    for d=1:2
        c_geod(d,:,k)=cumtrapz(linspace(0,1,T),v_geod(d,:,k));
    end
end

for k=1:S
    figure(2),clf;
    
    plot(c_geod(1,:,k),c_geod(2,:,k))
    title('Geodesic (click through steps)')
    axis equal    
    pause
end