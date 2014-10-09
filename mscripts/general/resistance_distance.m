

function R = resistance_distance(G)
%% Compute effective resistances with a matrix inversion            %%
%% Reference: D. Babi´c, D. J. Klein, I. Lukovits, S. Nikoli´c, and %%
%% N. Trinajsti´c. Resistance-distance matrix: A com-putational     %%
%% algorithm and its application. International Journal of Quantum  %%
%% Chemistry, 90(1):166–176, 2002                                   %%


% G is the matrix of conductances
% L is the graph Laplacian
n = nrow(G);
L = diag(sum(G))-G;
% H is the inverse of the sum-matrix
% H = inv(L + 1);
H = mldivide(L + 1,eye(n));
h1t = repmat(diag(H),1,n);
R = h1t + h1t' - 2*H;
