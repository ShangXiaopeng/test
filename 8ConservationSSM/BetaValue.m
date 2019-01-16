function [beta1,beta2,beta3,beta4] = BetaValue(x,NS)


x=x-6;

beta1=zeros(NS,NS,NS);
beta2=zeros(NS,NS);
beta3=zeros(NS,1);
beta4=zeros(NS,NS);



for L=1:NS
            
    XL1=x(L);
    XL2=x(L+1);
    V1=pi/6*10^(3*XL1);    
        
    if L>1
        
        for i1=1:L-1
            
            XLow=x(i1);
            XUp=x(i1+1);
            VI=pi/6*10^(3*XUp);
            
            for j1=1:L-1
                
                YLow=x(j1);
                YUp=x(j1+1);
                UJ=pi/6*10^(3*YUp);
                
%                 if (i1<L-1)&&(j1<L-1)
%                     
%                     beta1(i1,j1,L)=0;
%                     
%                 else
                    if XL1==YUp
                        a=XLow;
                    else
                        a=max(XLow,1/3*log10(10^(3*XL1)-10^(3*YUp)));
                    end
                    b=XUp;
                    
                    if (XL1==XUp)
                        
                        c=YLow;
                    else
                        c=max(YLow,1/3*(log10(10^(3*XL1)-10^(3*XUp))));
                    end
                    d=YUp;
                    
%                     beta1(i1,j1,L)=integral2(@ (X,Y) Beta1Fun(X,Y,XL1,XL2,XLow,XUp,YLow,YUp),XLow,XUp,YLow,YUp,'Method','iterated');

                    if (VI+UJ<V1)   beta1(i1,j1,L)=0;
                    else
                        beta1(i1,j1,L)=integral2(@ (X,Y) Beta1Fun(X,Y,XL1,XL2,XLow,XUp,YLow,YUp),a,b,c,d,'Method','iterated');
                    end
                    
%                 end
                
            end
                        
            a=XLow;
            b=XUp;
            c=max(XL1,1/3*(log10(10^(3*XL2)-10^(3*XUp))));
            d=XL2;
            
%             beta2(i1,L)=integral2(@ (X,Y) Beta2Fun(X,Y,XLow,XUp,XL1,XL2),XLow,XUp,XL1,XL2,'Method','iterated');
            beta2(i1,L)=integral2(@ (X,Y) Beta2Fun(X,Y,XLow,XUp,XL1,XL2),a,b,c,d,'Method','iterated');

            
        end
        
    end
    
    
    
    beta3(L)=integral2(@ (X,Y) Beta3Fun(X,Y,XL1,XL2),XL1,XL2,XL1,XL2,'Method','iterated');
%     beta3(L)=dblquad(@ (X,Y) Beta3FunHI(X,Y,XL1,XL2),XL1,XL2,XL1,XL2);
%     beta3(L)=Integ3 (m,n,XL1,XL2);
    
%     sum3=beta3;
    
    
    if L<NS
        
        for i2=L+1:NS
            
            XLow=x(i2);
            XUp=x(i2+1);
            
            beta4(i2,L)=integral2(@ (X,Y) Beta4Fun(X,Y,XL1,XL2,XLow,XUp),XLow,XUp,XL1,XL2,'Method','iterated');
%             beta4(i2,L)=dblquad(@ (X,Y) Beta4FunHI(X,Y,XL1,XL2,XLow,XUp),XLow,XUp,XL1,XL2);
%             beta4(i2,L)= Integ4 (m,n,XL1,XL2,XLow,XUp);
            
            
        end
        
    end




end
end