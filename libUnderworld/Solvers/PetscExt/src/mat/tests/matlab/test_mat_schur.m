%%
%%  test mat schur
%%

format long;


% ------------------------------------------------
%%  test_schur_11_test1
disp('***   test_schur_11_test1   ***');
clear;

M = 3;
N = 2;

K = MatIdentityNew(M);
G = MatFillStride( M,N, 2.0, 4.4 );
D = G';

p = ones(N,1);

disp('MatMult');
S = D * inv(K) * G;
fhat =  S * p
disp('MatMultTranspose');
St = S';
fhat = St * p


% ------------------------------------------------
%%  test_schur_11_test2
disp('***   test_schur_11_test2   ***');
clear;

M = 3;
N = 2;

K = MatIdentityNew(M);
G = MatFillStride( M,N, 2.0, 4.4 );
D = MatFillStride( N,M, 5.0, 1.3 );
C = MatFillStride( N,N, 1.0, 1.0 );

p = ones(N,1);


disp('MatMult');
S = D * inv(K) * G - C;
fhat = S * p
disp('MatMultTranspose');
St = S';
fhat = St * p


% ------------------------------------------------
%%  test_schur_22_test3 (& test_schur_22_test4)
disp('***   test_schur_22_test3  (& test_schur_22_test4) ***');
clear;

M = 3;
N = 2;

K = MatIdentityNew(M);
G = MatFillStride( M,N, 2.0, 4.4 );
D = MatFillStride( N,M, 5.0, 1.3 );
C = MatIdentityNew(N);
C = 2.0 * C;
alpha = 4.9;

u = ones(M,1);


disp('MatMult');
S = alpha * ( G * inv(C) * D - K );
fhat = S * u
disp('MatMultTranspose');
St = S';
fhat = St * u
