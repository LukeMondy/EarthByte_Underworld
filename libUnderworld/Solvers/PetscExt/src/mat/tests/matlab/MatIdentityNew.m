% ------------------------------------------------------------------
% I know there is a better way to do this, but it makes the matlab
% code and the c code look more similar if I have this function.
function A = MatIdentityNew( M )

A = zeros(M,M);
for i= 1:M
    A(i,i) = 1.0;
end
