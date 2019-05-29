function B=form_basis_path_straight(S,T,N1,N2,N3)

tgrid=linspace(0,1,T);
sgrid=linspace(0,1,S);

Tn = repmat(tgrid,S,1);
Sn = repmat(sgrid',1,T);

idx=1;

for s=1:N1
    for s1=1:N2
        B(:,:,idx) = ((cos(2*pi*s*Tn)-1).*(cos(2*pi*s1*Sn)-1))/(pi*sqrt(s*s1));
        B(:,:,idx+1) = ((cos(2*pi*s*Tn)-1).*sin(2*pi*s1*Sn))/(pi*sqrt(s*s1));
        B(:,:,idx+2) = (sin(2*pi*s*Tn).*sin(2*pi*s1*Sn))/(pi*sqrt(s*s1));
        B(:,:,idx+3) = (sin(2*pi*s*Tn).*(cos(2*pi*s1*Sn)-1))/(pi*sqrt(s*s1));
        idx=idx+4;
    end
end

for s1=1:N3
    B(:,:,idx) = (1-Tn).*sin(2*pi*s1*Sn)/(sqrt(pi)*sqrt(s1));
    B(:,:,idx+1) = (1-Tn).*(cos(2*pi*s1*Sn)-1)/(sqrt(pi)*sqrt(s1));
    B(:,:,idx+2) = Tn.*sin(2*pi*s1*Sn)/(sqrt(pi)*sqrt(s1));
    B(:,:,idx+3) = Tn.*(cos(2*pi*s1*Sn)-1)/(sqrt(pi)*sqrt(s1));

    idx=idx+4;
end

for j=1:size(B,3)
    B(end,:,j)=zeros(1,T);
end