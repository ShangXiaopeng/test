function beta4Fun = Beta4Fun(X,Y,XL1,XL2,XLow,XUp)

global kB T miu gama


v1=pi*10^(3*XL1)/6;
v2=pi*10^(3*XL2)/6;
u=pi*10.^(3*Y)/6;
v=pi*10.^(3*X)/6;



% v1=10.^XL1;
% v2=10.^XL2;
% u=10.^Y;
% v=10.^X;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

betaB=2*kB*T/3/miu*(u.^(1/3)+v.^(1/3)).*(u.^(-1/3)+v.^(-1/3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% beta4FunHI=beta_HY/1e12./v/(XUp-XLow)/(XL2-XL1);

beta=betaB;


beta4Fun=beta.*u.^gama; 



end

