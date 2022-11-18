function [xk, dif_norm, k, residuo] = Jacobi(A, b, x0, E, M, norma)
    
    U = triu(A, 1);
    L = tril(A, -1);
    D_inverse = diag(1./diag(A));
    
    xk = x0;
    k = 0;
    dif = E+1;
    
    while norm((dif), norma)>E & k<M
        xj = xk;
        xk = -(D_inverse)*(L+U)*xj + (D_inverse)*b;
        dif = xk-xj;
        k = k+1;
    end
    
    dif_norm = norm((dif), norma);
    residuo = norm((b-A*xk), norma);
    
endfunction
