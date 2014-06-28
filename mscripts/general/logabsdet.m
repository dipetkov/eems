function ld = logabsdet(A)
%% Compute log(abs(det(A))) %%


[L,U,P] = lu(A);
ld = sum(log(abs(diag(U))));
