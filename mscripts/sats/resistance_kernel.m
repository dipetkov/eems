

function kernel = resistance_kernel(Data,Graph,mValues,qValues)
%% Compute the auxiliary matrix X and the eigenvalues of X*inv(B) %%
%% necessary to compute efficiently the log likelihood terms      %%
%% Det(-inv(Delta)*Q) and tr(-inv(Delta)*Q*S)                     %%


G = sparse(Graph.Vi,Graph.Vj,mValues);
Ball = resistance_distance(G)/4;
Wall = reshape(qValues,[],1);

nSites = Data.nSites;
X = cell(nSites,1);
XC = cell(nSites,1);
ldDinv = zeros(nSites,1);
Bconst = zeros(nSites,1);
ldCiBwi = zeros(nSites,1);
cvtJtDJvct = cell(nSites,1);
oDinvoconst = zeros(nSites,1);

for s = 1:nSites
  Sizes = Data.Sizes{s};   %% the number of samples per deme
  oDemes = Data.oDemes{s}; %% the observed demes
  Cinv = diag(1./Sizes);
  C = diag(Sizes); 
  B = Ball(oDemes,oDemes);
  w = Wall(oDemes);
  winv = 1./w;
  BWinv = B*diag(winv);
  X{s} = mldivide(B*C-diag(w),BWinv);
  XC{s} = X{s}*C;
  ldCiBwi(s) = logabsdet(Cinv-BWinv);
  JtDJ = Data.JtDJ{s};
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

if Data.Testing
  oDemes = Data.oDemes2;
  Jinvpt = Data.Jinvpt2;
  n = length(Jinvpt);
  o = ones(n,1);
  Jpt = sparse(1:n,Jinvpt,1);
  B = Ball(oDemes,oDemes);
  w = Wall(oDemes);
  Delta = Jpt*B*Jpt' + Jpt*w*o'/2 + o*w'*Jpt'/2 - diag(Jpt*w);
  kernel.Delta = Delta;
end
