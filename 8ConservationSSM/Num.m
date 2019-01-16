
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
clc

set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'DefaultAxesFontName','Times New Roman') 



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

load('PHI1_1e-3.mat')

tau=t*beta0*N0;
phi_t=NT/N0;
avesize=1./phi_t;


Tau=0:0.0001:1;
Phi_t=PHI_S*(1+PHI_S*tanh(PHI_S*Tau/2))./(PHI_S+tanh(PHI_S*Tau/2));
AveSize=1./Phi_t;

plot(tau,avesize,'b-','LineWidth',1.5)
hold on
plot(Tau,AveSize,'b--','LineWidth',1.5)
hold on

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

load('PHI05_1e-3.mat')

tau=t*beta0*N0;
phi_t=NT/N0;
avesize=1./phi_t;


Tau=0:0.0001:1;
Phi_t=PHI_S*(1+PHI_S*tanh(PHI_S*Tau/2))./(PHI_S+tanh(PHI_S*Tau/2));
AveSize=1./Phi_t;

plot(tau,avesize,'r-','LineWidth',1.5)
hold on
plot(Tau,AveSize,'r--','LineWidth',1.5)
hold on

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

load('PHI10_1e-3.mat')


tau=t*beta0*N0;
phi_t=NT/N0;
avesize=1./phi_t;


Tau=0:0.0001:1;
Phi_t=PHI_S*(1+PHI_S*tanh(PHI_S*Tau/2))./(PHI_S+tanh(PHI_S*Tau/2));
AveSize=1./Phi_t;

plot(tau,avesize,'k-','LineWidth',1.5)
hold on
plot(Tau,AveSize,'k--','LineWidth',1.5)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 




xtick=0:0.2:1;
set(gca,'XTick',xtick,'XTickLabel',{'0', '0.2', '0.4', '0.6', '0.8', '1.0'})
ytick=0:0.2:1.4;
set(gca,'YTick',ytick,'YTickLabel',{'0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.2', '1.4'})


xlabel('\fontsize{11}\tau')
ylabel('\fontsize{11}1/\Phi(\tau)')



axis([-0.05 1.05 0 1.4])

set(gca,'FontSize',11,'LineWidth',1.0,'TickLength',[0.02 0.02])

set(gcf,'Units','centimeters','Position',[10 8 9.1 9.1]);
set(gca,'Position',[0.145 0.125 0.82 0.85])







% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



