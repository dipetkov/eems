

function kernel = resistance_kernel(Sstruct,Mstruct,mValues,qValues)
%% Compute the auxiliary matrix X and the eigenvalues of X*inv(B) %%
%% necessary to compute efficiently the log likelihood terms      %%
%% Det(-inv(Dhat)*Q) and tr(-inv(Dhat)*Q*S)                       %%


G = sparse(Mstruct.Mi,Mstruct.Mj,mValues);
Ball = resistance_distance(G)/4;
Wall = reshape(qValues,[],1);

nSites = Sstruct.nSites;
X = cell(nSites,1);
XC = cell(nSites,1);
ldDinv = zeros(nSites,1);
Bconst = zeros(nSites,1);
ldCiBwi = zeros(nSites,1);
oDinvoconst = zeros(nSites,1);

for s = 1:nSites
  Sizes = Sstruct.Sizes{s};   %% the number of samples per deme
  oDemes = Sstruct.oDemes{s}; %% the observed demes
  Cinv = diag(1./Sizes);
  C = diag(Sizes); 
  B = Ball(oDemes,oDemes);
  w = Wall(oDemes);
  winv = 1./w;
  BWinv = B*diag(winv);
  X{s} = mldivide(B*C-diag(w),BWinv);
  XC{s} = X{s}*C;
  ldCiBwi(s) = logabsdet(Cinv-BWinv);
  JtDJ = Sstruct.JtDJ{s};
  cwinvt = Sizes*winv';
  Bconst(s) = winv'*JtDJ*winv;            %% winv'*J'*D*J*winv
  oDinvoconst(s) = - Sizes'*winv;         %% -diag(J*J')'*winv
  ldDinvconst(s) = log(sum(Sizes)) ...
		 - sum(log(Sizes)) ...    %% logdet(C) = logdet(diag(Sizes))
		 + sum(log(winv).*Sizes); %% logdet(W) = logdet(diag(J*w))
  cvtJtDJvct{s} = JtDJ*cwinvt' + cwinvt*JtDJ';
end

kernel = struct('X',{X},'XC',{XC},...
		'Bconst',{Bconst},...
	 	'ldCiBwi',{ldCiBwi},...
		'cvtJtDJvct',...
		{cvtJtDJvct},...
		'ldDinvconst',...
		{ldDinvconst},...
		'oDinvoconst',...
		{oDinvoconst});
