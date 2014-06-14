function A = MatFillStride( M, N, start, stride )

% Return matrix
%   arg 1,2 - row col of matrix 
%   arg 3   - start_val
%   arg 4   - stride
%


A = zeros(M,N);

inc = 0.0;
for i= 1:M
    for j= 1:N
        A(i,j) = start + inc;
        inc = inc + stride;
    end
end

