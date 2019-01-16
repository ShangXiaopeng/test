function n= InitialFun(x)

global N0 V

v=10.^x;
% n0=N0*(2*N0/V)^2*v.*exp(-2*N0/V*v).*v*log(10);
n0=N0^2/V*exp(-N0/V*v).*v*log(10);


n=n0;

end





