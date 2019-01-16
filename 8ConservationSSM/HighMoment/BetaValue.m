function [beta1,beta2,beta3,beta4] = BetaValue(x,NS)

global gama

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
            VLow=pi*10^(3*XLow)/6;
            VUp=pi*10^(3*XUp)/6;
            
            for j1=1:L-1
                
                YLow=x(j1);
                YUp=x(j1+1);
                ULow=pi*10^(3*YLow)/6;
                UUp=pi*10^(3*YUp)/6;
                
                if (VUp+UUp<V1)   
                    
                    beta1(i1,j1,L)=0;
                    
                else
                    
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
                     beta1(i1,j1,L)=integral2(@ (X,Y) Beta1Fun(X,Y,XL1,XL2,XLow,XUp,YLow,YUp),a,b,c,d,'Method','iterated');
                     beta1(i1,j1,L)=beta1(i1,j1,L)/(VUp^gama-VLow^gama)/(UUp^gama-ULow^gama);    
                end
                    
%                 end
                
            end
                        
            a=XLow;
            b=XUp;
            c=max(XL1,1/3*(log10(10^(3*XL2)-10^(3*XUp))));
            d=XL2;
            
            ULow=pi*10^(3*XL1)/6;
            UUp=pi*10^(3*XL2)/6;
            VLow=pi*10^(3*XLow)/6;
            VUp=pi*10^(3*XUp)/6;
            
%             beta2(i1,L)=integral2(@ (X,Y) Beta2Fun(X,Y,XLow,XUp,XL1,XL2),XLow,XUp,XL1,XL2,'Method','iterated');
            beta2(i1,L)=integral2(@ (X,Y) Beta2Fun(X,Y,XLow,XUp,XL1,XL2),a,b,c,d,'Method','iterated');
            beta2(i1,L)=beta2(i1,L)/(VUp^gama-VLow^gama)/(UUp^gama-ULow^gama);

            
        end
        
    end
    
    
    ULow=pi*10^(3*XL1)/6;
    UUp=pi*10^(3*XL2)/6;
    VLow=pi*10^(3*XL1)/6;
    VUp=pi*10^(3*XL2)/6;
    
    beta3(L)=integral2(@ (X,Y) Beta3Fun(X,Y,XL1,XL2),XL1,XL2,XL1,XL2,'Method','iterated');
    beta3(L)=beta3(L)/(VUp^gama-VLow^gama)/(UUp^gama-ULow^gama);
    
    
    if L<NS
        
        for i2=L+1:NS
            
            XLow=x(i2);
            XUp=x(i2+1);
            ULow=pi*10^(3*XL1)/6;
            UUp=pi*10^(3*XL2)/6;
            VLow=pi*10^(3*XLow)/6;
            VUp=pi*10^(3*XUp)/6;

            beta4(i2,L)=integral2(@ (X,Y) Beta4Fun(X,Y,XL1,XL2,XLow,XUp),XLow,XUp,XL1,XL2,'Method','iterated');
            beta4(i2,L)=beta4(i2,L)/(VUp^gama-VLow^gama)/(UUp^gama-ULow^gama);
            
        end
        
    end




end
end