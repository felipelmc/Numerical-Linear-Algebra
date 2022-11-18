D = diag(fix(rand(10, 1, "uniform")*1000));
R = fix(rand(10, 10, "uniform")*10);
A = D + R;
x = fix(rand(10, 1, "uniform")*10);
b = A*x;
x0 = zeros(10, 1);

D = diag(fix(rand(100, 1, "uniform")*100000));
R = fix(rand(100, 100, "uniform")*100);
A = D + R;
x = fix(rand(100, 1, "uniform")*100);
b = A*x;
x0 = zeros(100, 1);

D = diag(fix(rand(1000, 1, "uniform")*100000));
R = fix(rand(1000, 1000, "uniform")*100);
A = D + R;
x = fix(rand(1000, 1, "uniform")*100);
b = A*x;
x0 = zeros(1000, 1);

D = diag(fix(rand(2000, 1, "uniform")*100000));
R = fix(rand(2000, 2000, "uniform")*100);
A = D + R;
x = fix(rand(2000, 1, "uniform")*100);
b = A*x;
x0 = zeros(2000, 1);

D = diag(fix(rand(5000, 1, "uniform")*100000));
R = fix(rand(5000, 5000, "uniform")*100);
A = D + R;
x = fix(rand(5000, 1, "uniform")*100);
b = A*x;
x0 = zeros(5000, 1);

D = diag(fix(rand(10000, 1, "uniform")*100000));
R = fix(rand(10000, 10000, "uniform")*100);
A = D + R;
x = fix(rand(10000, 1, "uniform")*100);
b = A*x;
x0 = zeros(10000, 1);

tic
GaussSeidel_inv(A, b, x0, 0.01, 1000, 1);
toc

tic
GaussSeidel(A, b, x0, 0.01, 1000, 1);
toc
