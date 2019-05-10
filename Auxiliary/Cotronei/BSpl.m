% Computes uniform Bsplines 
function y = BSpl(m,x) 
y=0;
if m==1
    if x>=0 & x<1
        y=1;
    else
        y=0;
    end
end

a=zeros(1,500);

if m>=2
    for k=1:m-1
        a(k)=0;
        x1=x-k+1;
        if x1>= 0 & x1<1
            a(k)=x1;
        end
        if x1>=1 & x1<2
            a(k)=2-x1;
        end
    end
    
    for p=1:m-2
        for q=1:m-1-p
            a(q)=((x-q+1)* a(q)+(p+q+1-x)*a(q+1))/(p+1);
        end
    end
    
    y=a(1);
end
    
end


