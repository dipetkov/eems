

function kernel = resistance_kernel(Sstruct,Mstruct,mValues)
%% Compute the auxiliary matrix X and the eigenvalues of X*inv(B) %%
%% necessary to compute efficiently the log likelihood terms      %%
%% Det(-inv(Delta)*Q) and tr(-inv(Delta)*Q*S)                     %%


G = sparse(Mstruct.Mi,Mstruct.Mj,mValues);
L = diag(sum(G)) - G;

Cconst = Sstruct.Cconst;
Cinv = Sstruct.Cinv;
C = Sstruct.Counts;
o = Sstruct.oDemes;
J = Mstruct.Juniq;
K = Mstruct.Kuniq;

if Mstruct.allobsrv
  Hinv = L + 1/o;
else
  Hi11 = L(J,J) + 1/o;
  Hi22 = L(K,K) + 1/o;
  Hi12 = L(J,K) + 1/o;
  H22Hi21 = linsolve(Hi22,Hi12');
  Hinv = Hi11 - Hi12*H22Hi21;
end

W0 = ones(o,1)*Cconst;
W = diag(W0);
Winv = diag(1./W0);

Binv = -2*Hinv;
BinvW = Binv*W;
X = mldivide(C-BinvW,Winv);
ldDinv = logabsdet(X*BinvW);
kernel = struct('X',{X},'XC',{X*C},...
		'ldDinv',{ldDinv});
