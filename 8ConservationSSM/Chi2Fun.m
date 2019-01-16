function chi2 = Chi2Fun(X,Y,XL1,XL2,XLow,XUp)


% u=pi*10.^(3*Y)/6;
% v=pi*10.^(3*X)/6;

u=10.^Y;
v=10.^X;

global s0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% a1=1.142;
% a2=0.558;
% a3=0.999;
% 
% k=1.381e-23;
% T=298;
% miu=17.9e-6;
% lamdg=6.86e-8;  %% Unit m
% 
% Kn1=lamdg./(r1);
% Kn2=lamdg./(r2);
% 
% C1=1+Kn1.*(a1+a2*exp(-a3./Kn1));
% C2=1+Kn2.*(a1+a2*exp(-a3./Kn2));
% 
% betaB=2*k*T/3/miu*(u.^(1/3)+v.^(1/3)).*(C1.*u.^(-1/3)+C2.*v.^(-1/3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


chi=2./v;
gama=s0*v;
% f=1./(3*log(10)*u);

f=1./u/log(10);

chi2=chi.*gama./f/(XUp-XLow);

    




end

