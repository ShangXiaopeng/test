clear all
clc


NS=10;


% dx=2.1/NS;
% x=-1.05:dx:1.05;

 
x=[
   -3.0000
   -2.5000
   -2.0000
   -1.5000
   -1.0000
   -0.5000
         0
    0.5000
    1.0000
    1.5000
    2.0000
];
% x=x+6;
% 
% 
% x=x-6;

% v=pi/6*10.^(3*x);
v=10.^x;
% vm=10.^xm;


beta1=zeros(NS,NS,NS);
beta2=zeros(NS,NS);
beta3=zeros(NS,1);
beta4=zeros(NS,NS);


for L=1:NS
            
    XL1=x(L);
    XL2=x(L+1);
        
        
    if L>1
        
        for i1=1:L-1
            
            XLow=x(i1);
            XUp=x(i1+1);
            
            for j1=1:L-1
                
                YLow=x(j1);
                YUp=x(j1+1);
                
                if (i1<L-1)&&(j1<L-1)
                    
                    beta1(i1,j1,L)=0;
                    
                elseif ( (i1==L-1) && (j1<L-1) )
                    
                    beta1(i1,j1,L)=(v(j1+1)+v(j1))*(v(j1+1)-v(j1))/2;
                    
                elseif ( (i1==L-1) && (j1==L-1) )
                    
                    beta1(i1,j1,L)=(v(L)-v(L-1))^2-(v(L)-2*v(L-1))^2/2;
                    
                elseif ((i1<L-1) && (j1==L-1))
                    
                    beta1(i1,j1,L)=(v(i1+1)+v(i1))*(v(i1+1)-v(i1))/2;
                
                end
                
                beta1(i1,j1,L)=beta1(i1,j1,L)/(v(i1+1)-v(i1))/(v(j1+1)-v(j1));
                
            end
                        
            
            beta2(i1,L)=(v(i1+1)+v(i1))*(v(i1+1)-v(i1))/2;
            beta2(i1,L)=beta2(i1,L)/(v(i1+1)-v(i1))/(v(L+1)-v(L));

            
        end
        
    end
    
    
    
    beta3(L)=2*(v(L+1)-v(L))^2-(v(L+1)-2*v(L))^2+(v(L+1)-2*v(L))^2/2;
    beta3(L)=beta3(L)/(v(L+1)-v(L))/(v(L+1)-v(L));
    
    if L<NS
        
        for i2=L+1:NS
            
            XLow=x(i2);
            XUp=x(i2+1);
            
            
            beta4(i2,L)=1;
%             beta4(i2,L)=1/(v(i2+1)-v(i2))/(v(L+1)-v(L));

            
            
        end
        
    end


    

end



for jj=1:NS
    
    xlswrite(num2str(jj),beta1(:,:,jj));
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

xlswrite('beta2',beta2);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

xlswrite('beta3',beta3);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  

xlswrite('beta4',beta4);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

xlswrite('X',x);

