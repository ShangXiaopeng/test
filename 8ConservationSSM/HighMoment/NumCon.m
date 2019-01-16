
function dQ=NumCon(t,Q,beta1,beta2,beta3,beta4,x,NS)


global gama

dQ=zeros(NS,1);

    
for L=1:NS
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     

    sum1=0;
    sum2=0;
    sum4=0;
        
    if L>1
       
        for i1=1:L-1
                        
            for j1=1:L-1
                 
                sum1=sum1+beta1(i1,j1,L)*Q(i1)*Q(j1);
                
            end
                        
            sum2=sum2+beta2(i1,L)*Q(i1);
            
        end
        
    end
    
    sum3=beta3(L);
      
    if L<NS
        
        for i2=L+1:NS
            
            sum4=sum4+beta4(i2,L)*Q(i2);
            
        end
        
    end




    dQ(L,1)=(3*gama*log(10))^2*(1/2*sum1-Q(L)*sum2-1/2*sum3*Q(L)^2-Q(L)*sum4);
   
    
end


disp(t);
end







