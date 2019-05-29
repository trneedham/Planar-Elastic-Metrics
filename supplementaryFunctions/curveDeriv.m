function v=curveDeriv(c)

[d,n]=size(c);

for i=1:d
    v(i,:)=gradient(c(i,:),1/(n));
end