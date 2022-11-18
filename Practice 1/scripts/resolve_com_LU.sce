function [X]=Resolve_com_LU(P, C, B) 
U = triu(C);
L = tril(C);
[n]=size(C,1);

for i=1:n
    L(i, i) = 1;
end

X = [];     

for i = 1 : n
    b = B(:, i);
    x = inv(U)*inv(L)*P*b;
    X = [X, x];
end
endfunction 
