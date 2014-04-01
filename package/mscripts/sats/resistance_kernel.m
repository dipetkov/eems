

function kernel = resistance_kernel(Sstruct,Mstruct,mValues)
%% Compute the auxiliary matrix X and the eigenvalues of X*inv(B) %%
%% necessary to compute efficiently the log likelihood terms      %%
%% Det(-inv(Delta)*Q) and tr(-inv(Delta)*Q*S)                     %%


Ball = sparse_solve_for_R(Mstruct,mValues);

nSites = Sstruct.nSites;
X = cell(nSites,1);
XC = cell(nSites,1);
ldDinv = zeros(nSites,1);

for s = 1:nSites

  C = Sstruct.Counts{s};
  Cinv = Sstruct.Cinv{s};
  Juniq = Mstruct.Juniq{s};
  B = Ball(Juniq,Juniq);
  o = Sstruct.oDemes(s);

  W = eye(o);
  BWinv = B;
  X{s} = mldivide(B*C-W,BWinv);
  XC{s} = X{s}*C;
  ldDinv(s) = - logabsdet(B-Cinv) ...
              - sum(log(diag(C)));

end

kernel = struct('X',{X},'XC',{XC},...
	 	'ldDinv',{ldDinv});
