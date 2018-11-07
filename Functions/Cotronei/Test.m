%Plot showing the growth of the polynomial reproduction coefficients
%associated to Daub or B-splines filters (but any scaling function can be
%chosen)

m=6; %sum-rule order (=order of poly reproduction) of the scaling function
% In the B-spline case, it corresponds to the order of the B-spline

a=-10; b=10; %Chosen iterval

% Daubechies
flag=0;
daub_filtername=['db' num2str(m)];
p=dbwavf(daub_filtername); p=p*2; %Daubechies filters of order m

%  B-spline of order m (uncomment if you wish to use it)
flag=1;
p=[ ];
for k=0:m
    p=[p;nchoosek(m,k)];
end
p=p*2^(-m+1);

N=length(p); %filter support [0,N-1]. In the B-spline case N-1 coincides with the order m
M=N-1;
close all


x=linspace(a,b,100);


for k=0:m-1
    % c=PolCoeffs(m,k,a); %coefficients that generate monomials of degree k
    c=PolCoeffs(m,k,p,a,b)
    
    % This applies only to B-splines
    if (flag==1)
        for i=1:length(x)
            pol(i)=0;
            jj=0;
            for l=-(M-a-1):b-1
                jj=jj+1;
                pol(i)=pol(i)+c(jj)*BSpl(m,x(i)-l);%
            end
        end
        figure(k+1)
        plot(x,x.^k,x,pol) %comparison between monomial of order k and the sum which reproduces it
        titolo2=['reproduction of the monomial of order ' num2str(k)'];
        title(titolo2)
    end
    
    figure(100)
    semilogy([-(M-a-1):b-1],abs(c),'*-')
    if k==0
        legenda=['k= ' num2str(k)];
        h=legend(legenda);
    end
    str=get(h,'string');
    legenda=  ['k= ' num2str(k)];
    h=legend([str legenda])  ;
    if flag==0
        scaling='Daubechies';
    else
        scaling='B-spline';
    end
    %  titolo=[ scaling ' of order M=' num2str(m) ' - growth of the polynomial reproduction coefficients'];
    %title(titolo)
    hold on
    
end

n=m-1;
for i=0:n
    cB=BernCoeffs(n,i,p,a,b);
    
    figure(300)
    plot([-(M-a-1):b-1],abs(cB),'*-')
    
    if i==0
        legenda=['i= ' num2str(i)];
        h=legend(legenda);
    end
    str=get(h,'string');
    legenda=  ['i= ' num2str(i)];
    h=legend([str legenda])  ;
    %      titolo=['Bernstein case - growth of the polynomial reproduction coefficients'];
    % title(titolo)
    hold on
    
    
    
    % This applies only to B-splines
    if (flag==1)
        
        for v=1:length(x)
            Bernst(v)=0;
            jj=0;
            for l=-(M-a-1):b-1
                jj=jj+1;
                Bernst(v)=Bernst(v)+cB(jj)*BSpl(m,x(v)-l); %
            end
        end
        figure(i+1+400)
        
        for t=1:length(x)
            B(t)=Bern1d(m-1,i,x(t),a,b);
        end
        
        plot(x,B,x,Bernst) %comparison between Bernstein poly of degree n and the sum which reproduces it
        legend('Exact','Reproduced')
        titolo2=['reproduction of the Bernstein polynomial of indexes n=' num2str(m-1) ' and i=' num2str(i)];
        title(titolo2)
    end
end

