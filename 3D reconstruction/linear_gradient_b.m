function [f, g] = linear_gradient_b(x,A,At,b)
% Generalized gradient step for ||Ax-b|| cost function
% Takes in variable, x, linear operator, A, and adjoint, At. Does not
% require passing in Atb
%
% outputs: 
% g: Gradient. Does not need to be a vector.
% f: Cost function value

u = A(x);
g = At(u-b);
f = norm(u(:)-b(:));