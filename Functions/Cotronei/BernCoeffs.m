function c=BernCoeffs(n,i,p,a,b)

% Computes the coefficients associated to the reproduction of the Bernstein polynomial  B_i^n
%  by the refinable function with mask p and with reproduction order m
% Interval: [a,b]
% Note that the filter p is normalized so that the sum of the filter coeffs is 2

M=length(p)-1;
l=0;

for h=-(M-a-1):b-1
    l=l+1;
    c(l)=0;
    
    if a==0
        for k=i:n
            pc=[];
            
            pc=PolCoeffs(n+1,k,p,a,b);
            
            c(l)=c(l)+nchoosek(n,k)*nchoosek(k,i)*(-1)^(k-i)*(1/b)^k*pc(l);
        end
    elseif b==0
        for k=n-i:n
            pc=[];
            
            pc=PolCoeffs(n+1,k,p,a,b);
            
            c(l)=c(l)+nchoosek(n,i)*nchoosek(i,k-n+i)*(-1)^(n-i)*(-1/a)^k*pc(l);
        end
        
    else
        
        for r=0:i
            for k=r:r+n-i
                pc=[];
                
                pc=PolCoeffs(n+1,k,p,a,b);
                
                c(l)=c(l)+nchoosek(i,r)*nchoosek(n-i,k-r)*(b/a)^(r)*(-1/b)^k*pc(l);
            end
        end
        c(l)=c(l)*nchoosek(n,i)*1/(b-a)^n*b^n*(-a/b)^i;
    end
    
end

