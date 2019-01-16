function mu= MuFun(x,Time)

global N0 DG0 SG0 gama kB T miu

dp=10.^(x);
v=pi/6*dp.^3*1e-18;

vg0=pi/6*DG0^3;
Z0=(log(SG0))^2;
K=2*kB*T/3/miu;
a=1+exp(Z0);
vg=vg0*(1+a*K*N0*Time)*exp(9/2*Z0)/(2+(exp(9*Z0)-2)/(1+a*K*N0*Time))^0.5;
dg=(6*vg/pi)^(1/3);
Z=1/9*log(2+(exp(9*Z0)-2)/(1+a*K*N0*Time));
Nt=N0/(1+a*K*N0*Time);

nt=log(10)*Nt/sqrt(2*pi)/Z^0.5.*exp(-(log(dp/dg)).^2/2/Z);

mu=v.^gama.*nt;


end





