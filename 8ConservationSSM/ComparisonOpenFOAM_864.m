
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all
clc
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')
set(0,'DefaultAxesFontName','Arial') 

% global miu kB T N0 DG0 SG0

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

zspan=[0,Time];

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


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % OpenFOAM Results Start% % % % % % % % % % % % % % % % 

QQ_OF=[
4.392460E-04
3.617030E-02
1.473150E+00
3.367420E+01
4.784540E+02
4.587790E+03
3.171910E+04
1.667520E+05
6.956120E+05
2.382610E+06
6.887170E+06
1.717570E+07
3.762180E+07
7.344730E+07
1.293640E+08
2.076840E+08
3.065530E+08
4.190050E+08
5.332830E+08
6.344430E+08
7.070120E+08
7.383090E+08
7.218080E+08
6.595130E+08
5.621220E+08
4.463460E+08
3.300850E+08
2.275610E+08
1.465280E+08
8.833610E+07
4.997830E+07
2.658810E+07
1.331770E+07
6.285310E+06
2.795880E+06
1.172260E+06
4.632350E+05
1.725020E+05
6.052670E+04
2.000850E+04
6.231100E+03
];

qq_OF=QQ_OF/dx;



% % % % % % % % % % % % % % % 

qq_OF(NS+1)=qq_OF(NS);
[XX,YY]=stairs(Dp,qq_OF);
figure
semilogx(XX,YY,'LineWidth',1.5)

hold on

% % % % % % % % % % % % % % % 

% B=bar(xm,qq_OF,0.7);
% B.FaceColor='b';
% 
% hold on

% % % % % % % % % % % % % % % 

% % % % % OpenFOAM Results End% % % % % % % % % % % % % % % % % 




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





















% semilogx(Dpm(1:NS),DQ0(1:NS),'r--d','LineWidth',2)
% hold on
% semilogx(Dpm(1:NS),q(1:NS),'b-+','LineWidth',2)
% hold on
% semilogx(Dpm(1:NS),q_Exp(1:NS)*1e12,'k-.+','LineWidth',2)



% set(gca,'LineWidth',1.0,'FontSize',12,'TickLength',[0.02 0.02])
% xlabel('\fontsize{14}Particle Diameter (\mum)')
% ylabel('\fontsize{14}\Delta\itN\rm/\Delta\itlog_{10}d\rm (number/m^{3})')
% xtick=10.^(-2:2);
% set(gca,'XTick',xtick,'XTickLabel',xtick)
% ytick=10.^(-1:0.1:0);
% ytick=0:2e11:1.2e12;

% set(gca,'YTick',ytick,'YTickLabel',{'0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'})
% set(gca,'YTick',ytick,'YTickLabel',{'0','2','4','6','8','10','12\times10^{11}'})
% axis([0.01,10,-0.5e11,1.2e12])
% legend('\fontsize{14}Initial PSD','\fontsize{14}Simulation','\fontsize{14}Experiment')


