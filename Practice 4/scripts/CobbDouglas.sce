function[alpha, b] = CobbDouglas(infos)
    bsys = log(infos(:, 1)) - log(infos(:, 3));
    A = [ones(bsys), log(infos(:, 2)) - log(infos(:, 3))];
    x = zeros(2);
    
    // aplicamos o Método dos Mínimos Quadrados
    // A_transpose * A * x = A_transpose * bsys
    x = Gaussian_Elimination_4(A'*A, A'*bsys);
    b = exp(x(1));
    alpha = x(2);
endfunction
    
    
    
    
    
