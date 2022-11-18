function [xk, dif_norm, k, residuo] = GaussSeidel_inv(A, b, x0, E, M, norma)
    
    U = triu(A, 1);
    L = tril(A);
    inverse = inv(L)
    
    xk = x0;
    k = 0;
    dif = E+1;
    
    while norm((dif), norma)>E & k<M
        xj = xk;
        xk = -inverse*U*xj + inverse*b;
        dif = xk-xj;
        k = k+1;
    end
    
    dif_norm = norm((dif), norma);
    residuo = norm((b-A*xk), norma);
    
endfunction
