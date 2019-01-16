
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
clc
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(0,'DefaultAxesFontName','Arial') 

global kB T miu 

N0=1e10;
DG0=1;
SG0=2.5;

P=101325;
kB=1.38054e-23;
T=293;
R=287;
gama=1.4;
miu0=1.71e-5;

% rhog=P/R/T;
miu=miu0*(273+111)/(T+111)*(T/273)^1.5;

Time=3600*24;


NS=41;
dx=4.2/NS;
x=-2.5:dx:1.7;

x=[
-2
-1.902439024
-1.804878049
-1.707317073
-1.609756098
-1.512195122
-1.414634146
-1.317073171
-1.219512195
-1.12195122
-1.024390244
-0.926829268
-0.829268293
-0.731707317
-0.634146341
-0.536585366
-0.43902439
-0.341463415
-0.243902439
-0.146341463
-0.048780488
0.048780488
0.146341463
0.243902439
0.341463415
0.43902439
0.536585366
0.634146341
0.731707317
0.829268293
0.926829268
1.024390244
1.12195122
1.219512195
1.317073171
1.414634146
1.512195122
1.609756098
1.707317073
1.804878049
1.902439024
2
];


xm=(x(1:NS)+x(2:NS+1))/2;
Dp=power(10,x);
Dpm=(Dp(1:NS)+Dp(2:NS+1))/2;


Q0=[
6227.900645
19992.85232
60454.21793
172186.0695
461944.4219
1167353.713
2778669.655
6230081.331
13157502.61
26174399.44
49046048
86567688.68
143923784.5
225390005.7
332477623.2
461972551
604639008.9
745422682.2
865636695.5
946881092.3
975623004.2
946881092.3
865636695.5
745422682.2
604639008.9
461972551
332477623.2
225390005.7
143923784.5
86567688.68
49046048
26174399.44
13157502.61
6230081.331
2778669.655
1167353.713
461944.4219
172186.0695
60454.21793
19992.85232
6227.900645    
];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

tspan=[0,Time];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

[beta1,beta2,beta3,beta4]=BetaValue(x,NS);

% [chi1,chi2,chi3]=ChiValue(x,NS);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

options=odeset('RelTol',1e-3,'MaxStep',100);

[t,Q]=ode45(@NumCon,tspan,Q0,options,beta1,beta2,beta3,beta4,x,NS);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

[M,N]=size(Q);
QQ=Q(M,:);
NT=sum(Q,2);



q=QQ./dx;
q(NS+1)=q(NS);

[XX,YY]=stairs(Dp,q);
figure
semilogx(XX,YY,'LineWidth',1.5)

hold on




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

DX=4.2/5/NS;
X=-2.5:DX:1.7;
DP=power(10,X);
% DPM=(DP(1:5*NS)+DP(2:5*NS+1))/2;
XM=(X(1:5*NS)+X(2:5*NS+1))/2;
DPM=power(10,XM);
vg0=pi/6*DG0^3;
Z0=(log(SG0))^2;
K=2*kB*T/3/miu;
a=1+exp(Z0);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % 

n0=log(10)*N0/sqrt(2*pi)/Z0^0.5.*exp(-(log(DPM/DG0)).^2/2/Z0);

semilogx(DPM,n0(1:5*NS),'k-','LineWidth',1.5)
% plot(log10(DPM),n0(1:5*NS),'k--','LineWidth',2.0)

hold on

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

vg=vg0*(1+a*K*N0*Time)*exp(9/2*Z0)/(2+(exp(9*Z0)-2)/(1+a*K*N0*Time))^0.5;
dg=(6*vg/pi)^(1/3);
Z=1/9*log(2+(exp(9*Z0)-2)/(1+a*K*N0*Time));
Nt=N0/(1+a*K*N0*Time);

nt=log(10)*Nt/sqrt(2*pi)/Z^0.5.*exp(-(log(DPM/dg)).^2/2/Z);

semilogx(DPM,nt(1:5*NS),'k--','LineWidth',1.5)
% plot(log10(DPM),nt(1:5*NS),'r--','LineWidth',2.0)


set(gca,'LineWidth',1.5,'fontname','arial','FontSize',9,'TickLength',[0.02 0.02])

xlabel('\fontname{arial}\fontsize{9}Particle Diameter (\mum)')
ylabel('\fontname{arial}\fontsize{9}\Delta\itN\rm/\Delta\itlog_{10}d_{p}\rm (number/m^{3}\times10^{-9})')

 set(gca,'XTickLabel',{'0.1','1.0','10'},'fontname','arial')
 set(gca,'YTickLabel',{'0','2','4','6','8','10','12'},'fontname','arial')
 
set(gcf,'Units','centimeters','Position',[10 8 10 8.5]);
set(gca,'Position',[0.15 0.15 0.8 0.8])

% axis([-1.1,1.25,0,1.2e10])

axis([0.05,20,0,1.2e10])



% legend('\fontsize{11}Simulation','\fontsize{11}Initial','\fontsize{11}Analytical')
% legend boxoff
% set(gca, 'FontName', 'Arial')



set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(0,'DefaultAxesFontName','Arial') 





