function[lambda1, x1, k, n_erro] = Potencia_deslocada_inversa(A, x0, epsilon, alfa, M)
    
    k=0;
    x0 = x0/norm(x0, 2);
    n_erro = epsilon + 1;
    
    NewA = A - (alfa*eye(A));
    
    while (k <= M && n_erro >= epsilon)
        x1 = Gaussian_Elimination_4(NewA, x0);
        x1 = x1/norm(x1, 2);
        lambda = x1'*A*x1;
        
        if (lambda < 0) then
            x1 = -x1;
        end
        
        n_erro = norm((x1-x0), 2);
        x0 = x1;
        k = k+1;
        
    end
   
   lambda1 = lambda
   
endfunction
