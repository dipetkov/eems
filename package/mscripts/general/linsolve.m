

function x = linsolve(A,b)
%% This is NOT the MATLAB linsolve function which takes options
%% to specify what kind of matrix A is
%% Here it is assumed that A is symmetric positive definite
%% and therefore - not triangular


%% Use a forward and a backward substitution 
%[R,p] = chol(A);
%x = R \ (R' \b);
x = mldivide(A,b);
