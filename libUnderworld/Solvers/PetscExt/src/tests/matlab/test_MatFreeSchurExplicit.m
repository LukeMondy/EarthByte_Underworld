%%
%%  test-MatFreeSchurExplicit
%%
clear;

format long;


%% schur_11
display( 'test_MatFreeSchurExplicit_11' );
M = 9;
N = 5;
G = MatFillStride(M,N, 1.0,2.2 );
D = MatFillStride(N,M, 2.0,1.2 );
C = MatFillStride(N,N, 1.0,0.0 );

% K is just diagonal so we leave it out

S11 = D*G-C


%% schur_22
display( 'test_MatFreeSchurExplicit_22' );
M = 9;
N = 5;
G = MatFillStride(M,N, 1.0,2.2 );
D = MatFillStride(N,M, 2.0,1.2 );
K = MatFillStride(M,M, 3.0,5.76 );

% C is just diagonal so we leave it out

S22 = G*D-K
