function [chi1,chi2,chi3] = ChiValue(x,NS)

% x=x-6;

chi1=zeros(NS,1);
chi2=zeros(NS,NS);
chi3=zeros(NS,1);

for L=1:NS
            
    XL1=x(L);
    XL2=x(L+1);
    
    chi1(L)=integral2(@ (X,Y) Chi1Fun(X,Y,XL1,XL2),XL1,XL2,XL1,XL2,'Method','iterated');
         
    if L<NS
        
        for i2=L+1:NS
            
            XLow=x(i2);
            XUp=x(i2+1);
            
            chi2(i2,L)=integral2(@ (X,Y) Chi2Fun(X,Y,XL1,XL2,XLow,XUp),XLow,XUp,XL1,XL2,'Method','iterated');
            
         end
        
    end

%     if L>1
        
%         chi3(L)=integral(@ (Y) Chi3Fun(Y,XL1,XL2),XL1,XL2,'Method','iterated');
        chi3(L)=integral(@ (Y) Chi3Fun(Y,XL1,XL2),XL1,XL2);


%     end

end
end