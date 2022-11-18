function [x, C, P]=Gaussian_Elimination_4(A, b)
    C=[A,b];
    [n]=size(C,1);
    P = eye(n,n);   // cria uma matriz identidade que servira para apontar as permutacoes
    for j=1:(n-1)
        maior_pivo = abs(C(j,j))       
        for i=j+1 : n
            if abs(C(i, j)) > maior_pivo then maior_pivo = C(i, j); // iteramos pelos elementos da matriz ate que seja encontrado um elemento maior que o pivo na posicao adequada
                a=i; // guarda na variavel "a" a linha do maior
            else a=j;
            end
        end
        C([j, a],:)= C([a, j],:);
        P([j, a],:)= P([a, j],:);  // realiza as mesmas permutacoes realizadas em A durante a eliminacao na matriz P, apontando quais permutacoes ocorreram


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
