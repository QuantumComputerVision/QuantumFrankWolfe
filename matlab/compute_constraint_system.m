function [A, b] = compute_constraint_system(n)

b = ones(n, 1);
I = eye(n);
A = [kron(I, b') ;  kron(b', I)] ;
b = [b; b];
end