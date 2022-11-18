function[lambda, x1, k, n_erro] = Metodo_potencia_1(A, x0, epsilon, M)
    
    k = 0;
    [valor, ind] = max(abs(x0));
    x0 = x0/x0(ind);
    x1 = A*x0;
    n_erro = epsilon + 1;
    
    while (k <= M && n_erro >= epsilon)
        [valor, ind] = max(abs(x1));
        lambda = x1(ind);
        x1 = x1/lambda;
        n_erro = norm((x1-x0), %inf);
        x0 = x1;
        x1 = A*x0;
        k = k + 1;
    end
    
endfunction


function[lambda, x1, k, n_erro] = Metodo_potencia_2(A, x0, epsilon, M)
    
    k = 0;
    x0 = x0/norm(x0, 2);
    x1 = A*x0;
    n_erro = epsilon + 1;
    
    while (k <= M && n_erro >= epsilon)
        lambda = x1'*x0;
        
        if (lambda < 0) then
            x1 = -x1;
        end
        
        x1 = x1/norm(x1, 2);
        n_erro = norm((x1-x0), 2);
        x0 = x1;
        x1 = A*x0;
        k = k+1;
    end
    
endfunction
