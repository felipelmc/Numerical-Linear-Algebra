function[Q, R]= qr_GS(A)
    [m, n] = size(A); // le as dimensoes de A
    Q = zeros(m, n); // inicializa Q e R com zeros e as dimensoes de A
    R = zeros(n, n);
    
    for j = 1:n
    
        v = A(:, j); // v = a_j
        
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j); // r_ij = q_i^T * a_j
            v = v - R(i, j) * Q(:, i); // v = v - v_ij * q_i
        end
    
        R(j, j) = norm(v, 2); // r_jj = norma 2 de v
        Q(:, j) = v/R(j, j); // q_j = v/r_jj
    end
endfunction


function[Q, R]= qr_GSM(A)
    [m, n] = size(A) // le as dimensoes de A
    Q = zeros(m, n); // inicializa Q e R com zeros e as dimensoes de A
    R = zeros(n, n);
    
    for j = 1:n
        v = A(:, j); // v = a_j
        
        for i = 1:j-1
            R(i, j) = Q(:, i)' * v; // r_ij = q_i^T * v
            v = v - R(i, j) * Q(:, i); // v = v - v_ij * q_i
        end
        
        R(j, j) = norm(v, 2); // r_jj = norma 2 de v
        Q(:, j) = v/R(j, j); // q_j = v/r_jj
    end
endfunction


function[U, R] = qr_House(A) 
    [m, n] = size(A) // le as dimensoes de A
    k = min(m-1, n) // obtem o numero k de vetores necessarios 
    U = zeros(m, k) // inicializa U com zeros
    
    for j = 1:k
        x = A(j:m, j);
        
        if x(1) > 0
            x(1) = x(1) + norm(x, 2)        // calcula x - Hx
        else 
            x(1) = x(1) - norm(x, 2)
        end
        
        u = x / norm(x, 2)  // u = (x-Hx) / norma 2 de (x = Hx)
        U(j:m, j) = u       // guarda o vetor u em U
        A(j:m, j:n) = A(j:m, j:n) - 2 * u * (u' * A(j:m, j:n))
    end
    
    R = triu(A);
endfunction


function[Q] = constroi_Q_House(U)
    [m, k] = size(U)      // le as dimensoes de U
    Q = eye(m, m)       // inicializa Q como uma identidade mxm
    
    for j = 1:k
        u = U(1:m, j)
        Q = Q - 2 * (Q * u) * u'
    end
endfunction

function[S] = espectro(A, tol)
    [Q, R] = qr_GS(A) // decompoe A
    
    current = 0
    
    while current < tol // itera ate alcancar a tolerancia definida
        d1 = diag(A)
        A = Q' * A * Q
        d2 = diag(A)
        [Q, R] = qr_GS(A)
        current = norm(d1 - d2, %inf)
    end
    
    S = diag(A)
endfunction
