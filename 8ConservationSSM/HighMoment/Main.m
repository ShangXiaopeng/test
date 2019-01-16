
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
clc
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(0,'DefaultAxesFontName','Arial') 

global kB T miu gama N0 DG0 SG0

N0=1e10;
DG0=1;
SG0=2.5;

P=101325;
kB=1.38054e-23;
T=293;
gama=2/3;
miu0=1.71e-5;
miu=miu0*(273+111)/(T+111)*(T/273)^1.5;

Time=3600*24;

XL=-2;
XU=2.0;

NS=50;
dx=(XU-XL)/NS;
x=XL:dx:XU;

xm=(x(1:NS)+x(2:NS+1))/2;
Dp=power(10,x);
Dpm=(Dp(1:NS)+Dp(2:NS+1))/2;

v=pi*10.^(3*x)/6*1e-18;
vm=pi*10.^(3*xm)/6*1e-18;

for j=1:NS
    
    Q0(j)=integral(@InitialFun,x(j),x(j+1));
    
end
Q0=Q0./dx.*(v(2:NS+1).^gama-v(1:NS).^gama)/3/gama/log(10);

q0=Q0./dx;
q0(NS+1)=q0(NS);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

tspan=[0,Time];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

[beta1,beta2,beta3,beta4]=BetaValue(x,NS);

% [chi1,chi2,chi3]=ChiValue(x,NS);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

options=odeset('RelTol',1e-4,'MaxStep',100);

[t,Q]=ode45(@NumCon,tspan,Q0,options,beta1,beta2,beta3,beta4,x,NS);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

[M,N]=size(Q);
QQ=Q(M,:);
Mu_Sim=sum(Q,2);

q=QQ./(v(2:NS+1).^gama-v(1:NS).^gama);


% q(NS+1)=q(NS);
% 
% [XX,YY]=stairs(Dp,q);
% 
% semilogx(XX,YY)


plot(t,Mu_Sim,'b-','LineWidth',1.5)

hold on


NTime=10000;
TIME=0:Time/NTime:Time;

for j3=1:length(TIME)
    
    Mu(j3)=integral(@ (x) MuFun(x,TIME(j3)),XL,XU);
    
end

plot(TIME,Mu,'k--','LineWidth',1.5)



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% DX=4.0/5/NS;
% X=-2.0:DX:2.0;
% DP=power(10,X);
% % DPM=(DP(1:5*NS)+DP(2:5*NS+1))/2;
% XM=(X(1:5*NS)+X(2:5*NS+1))/2;
% DPM=power(10,XM);
% 
% vg0=pi/6*DG0^3;
% Z0=(log(SG0))^2;
% K=2*kB*T/3/miu;
% a=1+exp(Z0);
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% n0=log(10)*N0/sqrt(2*pi)/Z0^0.5.*exp(-(log(DPM/DG0)).^2/2/Z0);
% m0=(pi/6*DPM.^3*1e-18).^gama.*n0;
% 
% semilogx(DPM,n0(1:5*NS),'k-','LineWidth',1.5)
% 
% hold on
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% vg=vg0*(1+a*K*N0*Time)*exp(9/2*Z0)/(2+(exp(9*Z0)-2)/(1+a*K*N0*Time))^0.5;
% dg=(6*vg/pi)^(1/3);
% Z=1/9*log(2+(exp(9*Z0)-2)/(1+a*K*N0*Time));
% Nt=N0/(1+a*K*N0*Time);
% 
% nt=log(10)*Nt/sqrt(2*pi)/Z^0.5.*exp(-(log(DPM/dg)).^2/2/Z);
% mt=(pi/6*DPM.^3*1e-18).^gama.*nt;
% semilogx(DPM,nt(1:5*NS),'k--','LineWidth',1.5)
% 
% 
% 
% 
% 
% xlabel('\fontname{arial}\fontsize{9}Particle Diameter (\mum)')
% ylabel('\fontname{arial}\fontsize{9}\Delta\itN\rm/\Delta\itlog_{10}d_{p}\rm (number/m^{3}\times10^{-9})')
% 
% 
% set(gca,'LineWidth',1.5,'fontname','arial','FontSize',9,'TickLength',[0.02 0.02])
% set(gca,'XTickLabel',{'0.1','1.0','10'},'fontname','arial')
% set(gca,'YTickLabel',{'0','2','4','6','8','10','12'},'fontname','arial')
% set(gcf,'Units','centimeters','Position',[10 8 10 8.5]);
% set(gca,'Position',[0.15 0.15 0.8 0.8])



% axis([-1.1,1.25,0,1.2e10])
% axis([0.05,20,0,1.2e10])



% legend('\fontsize{11}Simulation','\fontsize{11}Initial','\fontsize{11}Analytical')
% legend boxoff
% set(gca, 'FontName', 'Arial')



set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(0,'DefaultAxesFontName','Arial') 





