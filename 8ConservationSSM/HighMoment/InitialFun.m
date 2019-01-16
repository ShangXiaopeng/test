function n= InitialFun(x)

global N0 DG0 SG0


% n0=N0*(2*N0/V)^2*v.*exp(-2*N0/V*v).*v*log(10);
% n0=N0^2/V*exp(-N0/V*v).*v*log(10);

n0=N0/sqrt(2*pi)/log(SG0)*log(10)*exp(-((log(10.^x/DG0)).^2/2/(log(SG0))^2));

n=n0;

end





