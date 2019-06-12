function [f]=ScalingIntegers(c)
%Computes  the scaling function at integers from the mask c

N=length(c)-1; %supporto=[0,N]
c=2*c/sum(c); %normalize to sum 2
for i=2:N
    for j=2:N
        if 2*i-j>=1 & 2*i-j<=N+1
            A(i-1,j-1)=c(2*i-j);
        end
    end
end
[V,D] = eig(A);

[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs = V(:,ind);

f=Vs(:,end);
f=[0;f;0]/sum(f); %valori normalizzati a 1 (partition of unity)

    


