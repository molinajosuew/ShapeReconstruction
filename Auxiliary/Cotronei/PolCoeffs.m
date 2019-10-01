function c=PolCoeffs(m,k,p,a,b)

% Computes the coefficients associated to the reproduction of monomials of
% degree k (k=0,...,m-1), in the interval [a,b], by the refinable function with mask p and with reproduction order m
% (i.e. satisfying sum-rules of order m)
% Ref. Chui, De Villiers, Wavelets subdivision methods, CRC Press

% Note that the filter p is normalized so that the sum of the filter coeffs is 2

M=length(p)-1;
qcoeffs=Qpol(m,p,m-1-k);
l=0;
for j=-(M-a-(mod(m,2)==0)):b-1 %The exact number of coefficients involved in the formula is b-a+M-1
    l=l+1;
%     c(l)=factorial(k)/factorial(m-1)*polyval(qcoeffs,j+(mod(m,2)==1)*.5);
    c(l)=factorial(k)/factorial(m-1)*polyval(qcoeffs,j);
    %c(j+M+1)=factorial(k)/factorial(m-1)*polyval(qcoeffs,j);
end
   
function [qcoeff]=Qpol(m,p,order)
%Computation of the polynomial Q (and of its derivatives up to order) as in Chui-De Villiers page 176  
mu=DiscreteMoments(m,p);

q=ones(m,1);
 for j=m-2:-1:0
     q(j+1)=0;
     for k=j+1:m-1
     q(j+1)=q(j+1)+(-1)^k*nchoosek(k,j)*mu(k+1-j)*q(k+1);
     end
     q(j+1)=q(j+1)*(-1)^(j+1);
 end
 
 q=q(end:-1:1);
 
 
 if(order==0)
     qder=q;
 else
      for k=1:order  
     qder=polyder(q);
     q=qder;
      end
 end

       qcoeff=qder(1:m-order);
     
function [mu]=DiscreteMoments(m,p)

       for k=0:m-1
    alpha(k+1)=0;
    for j=0:m
        if 2*j+1 <=length(p)
            alpha(k+1)=alpha(k+1)+(2*j)^k*p(2*j+1);
        end
    end
end

for k=0:m-1
    alpha2(k+1)=0;
    for j=0:m
        if 2*j+2  <=length(p)
            alpha2(k+1)=alpha2(k+1)+(2*j+1)^k*p(2*j+2);
        end
    end
end

if (abs(alpha-alpha2)>10^-7)
    error('Sum rules are not satisfied!');
end

mu=ones(m,1);
for k=1:m-1
    mu(k+1)=0;
    for j=0:k-1
        mu(k+1)=mu(k+1)+ nchoosek(k,j)*alpha(k+1-j)*mu(j+1);
    end
    mu(k+1)=mu(k+1)*1/(2^k-1);
end
