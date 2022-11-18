function [xk, dif_norm, k, residuo]=GaussSeidel(A, b, x0, E, M, norma) 
    
    U = triu(A, 1); 
    L = tril(A);

    xk = x0;
    xk_size = size(xk, 1); 
    k = 0; 
    dif = E + 1;

    while k<M & norm(dif, norma)> E
        x = xk; 
        U_b = -U*xk + b; 
        xk(1) = U_b(1)/L(1, 1);
        
        for i = 2:xk_size 
            xk(i) = (U_b(i) - L(i, 1:i-1)*xk(1:i-1))/L(i, i);
        end
    
        dif = xk - x 
        k = k+1;
    end

    dif_norm = norm((dif), norma);
    residuo = norm((b-A*xk), norma);
    
endfunction
