% Pseudo-spline filters and computation of polynomial reproduction coefficients.

close all
clear all

a=-2; b=2; %Chosen iterval
type=1; % Specify the type (0 for primals, 1 for duals. PLEASE NOTE THAT DUAL HERE DOES NOT MEAN DUAL
% IN THE SENSE OF WAVELETS. It is a notation used by
% people working in subdivision to specify that the center
% of the support is 1/2 instead of 0

m=8; %Specify  the order of poly reproduction (i.e. degree of poly repr. =m-1)
% m=2,4,6,8,10 for the dual pseudosplines
%  m=2,4,6,8 for the primal pseudosplines

if type==0
    switch m
        case 2
            p =[ 1/256, 9/256, 9/64, 21/64, 63/128, 63/128, 21/64, 9/64, 9/256, 1/256];
        case 4
            p = [ -9/2048, -55/2048, -99/2048, 99/2048, 363/1024, 693/1024, 693/1024, 363/1024, 99/2048, -99/2048, -55/2048, -9/2048];
        case 6
            p = [ 99/32768, 351/32768, -143/16384, -1287/16384, -1287/32768, ...
                10725/32768, 6435/8192, 6435/8192, 10725/32768, -1287/32768, ...
                -1287/16384, -143/16384, 351/32768, 99/32768];
        case 8
            p=[ -429/262144, -495/262144, 4095/262144, 5005/262144, -19305/262144, ...
                -27027/262144, 75075/262144, 225225/262144, 225225/262144, 75075/262144, ...
                -27027/262144, -19305/262144, 5005/262144, 4095/262144, -495/262144, ...
                -429/262144];
        case  10
            p=[ 6435/8388608, -7293/8388608, -8415/1048576, 9945/1048576, ...
                85085/2097152, -109395/2097152, -153153/1048576, 255255/1048576, ...
                3828825/4194304, 3828825/4194304, 255255/1048576, -153153/1048576, ...
                -109395/2097152, 85085/2097152, 9945/1048576, -8415/1048576, ...
                -7293/8388608, 6435/8388608];
            
    end
else
    switch m
        case 2
            p=1/128*[0, 1, 8, 28, 56, 70, 56, 28, 8, 1,0];
        case 4
            p=1/128*[0, -1, -5, -5, 20, 70, 98, 70, 20, -5, -5, -1,0 ];
        case 6
            p=1/1024*[0, 5, 12, -30, -100, 75, 600, 924, 600, 75, -100, -30, 12, 5,0];
        case 8
            p=1/4096*[0, -10, 0, 98, 0, -490, 0, 2450, 4096, 2450, 0, -490,0, 98,0, -10,0 ];
    end
end

N=length(p); %filter support [0,N-1].

offset=(N-1)/2; %IMPORTANT: this is the shift tau needed for evaluating correctly the coefficients as values of the monomials

p=2*p/sum(p); % Normalize so that sum=2

M=N-1;

%x=linspace(a,b,100);
z=[a:b];
[f]=ScalingIntegers(p);

% Computation and representation of monomials coefficients
for k=0:m-1
    
    [c,offset]=PseudosplinesPolCoeffs(m,k,p,a,b)
   
    for v=1:length(z)
        pol(v)=0;
        jj=0;
       
        for l=-(M-a-1):b-1
            jj=jj+1;
            
            if z(v)-l+1>=1 & z(v)-l+1<=length(f)
                pol(v)=pol(v)+c(jj)*f(z(v)-l+1)
                
            end
        end
        
    end
       
    figure(2220+k)
    
    plot([-(M-a-1):b-1],c,'o',[-(M-a-1):0.01:b-1],[-(M-a-1)+offset:0.01:b-1+offset].^k,'k')
    legend(['j^' num2str(k)],['x^' num2str(k)])
    title(['grado' num2str(k)])
    
    figure(k+1)
    plot(z,z.^k,z,pol) %comparison between monomial of order k and the sum which reproduces it
    titolo2=['reproduction of the monomial of degree ' num2str(k)'];
    legend('exact','reproduced')
    title(titolo2)
    
    
    figure(100)
    semilogy([-(M-a-1):b-1],abs(c),'*-')
    if k==0
        legenda=['k= ' num2str(k)];
        h=legend(legenda);
    end
    str=get(h,'string');
    legenda=  ['k= ' num2str(k)];
    h=legend([str legenda])  ;
    titolo=[ 'Pseudospline  of order m=' num2str(m) ' - growth of the polynomial reproduction coefficients'];
    title(titolo)
    hold on
    
end

% Computation and representation of Bernstein coefficients
n=m-1
for i=0:n
    cB=PseudosplinesBernCoeffs(n,i,p,a,b);

    figure(300)
    plot([-(M-a-1):b-1],abs(cB),'*-')

    if i==0
        legenda1=['i= ' num2str(i)];
        h1=legend(legenda1);
    end
    str=get(h1,'string');
    legenda1=  ['i= ' num2str(i)];
    h1=legend([str legenda1])  ;
    titolo=['Bernstein case - growth of the polynomial reproduction coefficients'];
    title(titolo)
    hold on
      
        for v=1:length(z)
            Bernst(v)=0;
            jj=0;
            for l=-(M-a-1):b-1
                jj=jj+1;
                if z(v)-l+1>=1 & z(v)-l+1<=length(f)
                    Bernst(v)=Bernst(v)+cB(jj)*f(z(v)-l+1);%
                end
            end
        end
        figure(i+1+400)

        for t=1:length(z)
            B(t)=Bern1d(m-1,i,z(t),a,b);
        end


        plot(z,B,z,Bernst) %comparison between Bernstein poly of degree n and the sum which reproduces it
        legend('Exact','Reproduced')
        titolo2=['reproduction of the Bernstein polynomial of indexes n=' num2str(m-1) ' and i=' num2str(i)];
        title(titolo2)
     end
