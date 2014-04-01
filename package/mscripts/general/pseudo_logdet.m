function ld = pseudo_logdet(A,r)
%% log(Det(A)) where Det is the pseudo determinant %%
%% [ie, the product of all non-zero eigenvalues]   %%


tol = 1e-12;
if nargin==1
  r = rank(A);
end

D = eig(A);
D = sort(D,'descend');
D = D(1:r);
ld = sum(log(D));

if abs(imag(ld))<tol
  ld = real(ld);
end
