function [x, C]=Gaussian_Elimination_3(A, b)
C=[A,b];
[n]=size(C,1);
for j=1:(n-1)
    if C(j,j) == 0
    then
    maior_pivo = abs(C(j,j))       
    for i=j+1 : n
        if abs(C(i, j)) > maior_pivo then maior_pivo = C(i, j);
          a=i;
        else a=j;
        end
    end
    C([j, a],:)= C([a, j],:);
end

for i=(j+1):n
 C(i,j)=C(i,j)/C(j,j);
 C(i,j+1:n+1)=C(i,j+1:n+1)-C(i,j)*C(j,j+1:n+1);;
end
end

x=zeros(n,1);
x(n)=C(n,n+1)/C(n,n);
for i=n-1:-1:1
 x(i)=(C(i,n+1)-C(i,i:n)*x(i:n))/C(i,i);
end
C=C(1:n,1:n);
endfunction
