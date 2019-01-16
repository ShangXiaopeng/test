
clear all
clc

global N0 DG0 SG0

N0=1e10;
DG0=1;
SG0=2.5;

NS=41;
dx=4/NS;

x=-2.0:dx:2.0;

Dp=power(10,x);

Dpm=(Dp(1:NS)+Dp(2:NS+1))/2;

% n=N0/(2*pi)^0.5/log(SG)./Dp.*exp(-(log(Dp)-log(DG)).^2/2/log(SG)^2);
% loglog(Dp,n)


for j=1:NS
    
    Q(j)=quad(@InitialFun,x(j),x(j+1));
    
end

Q(NS+1)=Q(NS);
q=Q./dx;



[XX,YY]=stairs(Dp,q);

semilogx(XX,YY)

hold on

semilogx(Dpm,q(1:NS))

grid on

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

xlswrite('X',x');

xlswrite('Q',Q(1:NS)');




