
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

NSS=1000;
X=x1:(x2-x1)/NSS:x2;
XM=(X(1:NSS)+X(2:NSS+1))/2;
VM=10.^XM;

Eta=N0*VM/V;

PHI=PHI_S*(1+PHI_S*tanh(PHI_S*TAU/2))/(PHI_S+tanh(PHI_S*TAU/2));

Phi0=exp(-Eta).*Eta*log(10);
% Phi0=exp(-Eta);


Phi=PHI^2*exp(-Eta.*PHI).*Eta*log(10);
% Phi=PHI^2*exp(-Eta.*PHI);

semilogx(Eta,Phi0,'b--','LineWidth',2)

hold on

semilogx(Eta,Phi/PHI^2,'b-','LineWidth',2)



% hold on
% phi0=Q0./dv*V/N0^2.*eta*log(10);
% semilogx(eta,phi0,'bs')



hold on
phi_t=sum(QQ)/sum(Q0);
phi=QQ./dv*V/N0^2.*eta*log(10);
semilogx(eta,phi/phi_t^2,'bs')
hold on

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

load('PHI05_1e-3.mat')

NSS=1000;
X=x1:(x2-x1)/NSS:x2;
XM=(X(1:NSS)+X(2:NSS+1))/2;
VM=10.^XM;

Eta=N0*VM/V;

PHI=PHI_S*(1+PHI_S*tanh(PHI_S*TAU/2))/(PHI_S+tanh(PHI_S*TAU/2));

Phi0=exp(-Eta).*Eta*log(10);
% Phi0=exp(-Eta);


Phi=PHI^2*exp(-Eta.*PHI).*Eta*log(10);
% Phi=PHI^2*exp(-Eta.*PHI);

semilogx(Eta,Phi0,'r--','LineWidth',2)

hold on

semilogx(Eta,Phi/PHI^2,'r-','LineWidth',2)


% hold on
% phi0=Q0./dv*V/N0^2.*eta*log(10);
% semilogx(eta,phi0,'ro')



hold on
phi_t=sum(QQ)/sum(Q0);
phi=QQ./dv*V/N0^2.*eta*log(10);
semilogx(eta,phi/phi_t^2,'ro')
hold on

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


load('PHI10_1e-3.mat')

NSS=1000;
X=x1:(x2-x1)/NSS:x2;
XM=(X(1:NSS)+X(2:NSS+1))/2;
VM=10.^XM;

Eta=N0*VM/V;

PHI=PHI_S*(1+PHI_S*tanh(PHI_S*TAU/2))/(PHI_S+tanh(PHI_S*TAU/2));

Phi0=exp(-Eta).*Eta*log(10);
% Phi0=exp(-Eta);


Phi=PHI^2*exp(-Eta.*PHI).*Eta*log(10);
% Phi=PHI^2*exp(-Eta.*PHI);

semilogx(Eta,Phi0,'k--','LineWidth',2)

hold on

semilogx(Eta,Phi/PHI^2,'k-','LineWidth',2)



% hold on
% phi0=Q0./dv*V/N0^2.*eta*log(10);
% semilogx(eta,phi0,'k^')



hold on
phi_t=sum(QQ)/sum(Q0);
phi=QQ./dv*V/N0^2.*eta*log(10);
semilogx(eta,phi/phi_t^2,'k^')





% xtick=0:0.2:1;
% set(gca,'XTick',xtick,'XTickLabel',{'0', '0.2', '0.4', '0.6', '0.8', '1.0'})
xtick=10.^(-3:2);
set(gca,'XTick',xtick,'XTickLabel',{'0.001', '0.01', '0.1', '1.0', '10', '100'})
ytick=0:0.2:1.2;
set(gca,'YTick',ytick,'YTickLabel',{'0', '0.2', '0.4', '0.6', '0.8', '1.0', '1.2'})


xlabel('\fontsize{13}\eta')
ylabel('\fontsize{13}\phi(\eta,\tau=1)')



axis([0.0009 110 -0.05 1.205])


set(gca,'FontSize',13,'LineWidth',1.0,'TickLength',[0.02 0.02])

set(gcf,'Units','centimeters','Position',[10 8 12.1 9.0]);
set(gca,'Position',[0.145 0.125 0.82 0.85])




