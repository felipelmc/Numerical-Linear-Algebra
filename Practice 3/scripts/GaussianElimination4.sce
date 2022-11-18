function [x, C, P]=Gaussian_Elimination_4(A, b)

    C=[A,b];
    [n]=size(C,1);
    P = eye(n,n);

    for j=1:(n-1)

        [valor, ind] = max(abs(C([j:n], j))); 

        if j ~= ind then
            c_ind = ind + j - 1;
            C([j, c_ind],:)= C([c_ind, j],:);
            P([j, c_ind],:)= P([c_ind, j],:);
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
