% Takes a curve p and preprocesses it by: 
% -shifting to begin at the origin
% -scaling so that it has total length 1.

function pnew=preprocess_curve(p)

[d,n]=size(p);

% Shift to begin at the origin.
for i=1:d
   pnew2(i,:)=p(i,:)-p(i,1);
end

pnew2d = gradient(pnew2,1/n);
 
curveLength = trapz(linspace(0,1,n),sqrt(sum(pnew2d.^2,1)));

% Rescale so that the curve has total length 1.

pnew=pnew2/curveLength;