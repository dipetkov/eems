

function kernel = resistance_kernel(Data,Graph,mValues,qValues)
%% Compute the auxiliary matrix X and the eigenvalues of X*inv(B) %%
%% necessary to compute efficiently the log likelihood terms      %%
%% Det(-inv(Delta)*Q) and tr(-inv(Delta)*Q*S)                     %%


G = sparse(Graph.Vi,Graph.Vj,mValues);
Ball = resistance_distance(G);
Wall = reshape(qValues,[],1);

nSites = Data.nSites;
X = cell(nSites,1);
Xc_winv = cell(nSites,1);
oDinvo = zeros(nSites,1);
ldCiBwi = zeros(nSites,1);

for s = 1:nSites
  Sizes = Data.Sizes{s};   %% the number of samples per deme
  oDemes = Data.oDemes{s}; %% the observed demes
  B = Ball(oDemes,oDemes);
  w = Wall(oDemes)*2;
  winv = 1./w;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [L,U,P] = lu(B*diag(Sizes) - diag(w));
  X{s} = mldivide(U,mldivide(L,P*B*diag(winv)));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Bc = B*diag(Sizes);
  [L,U,P] = lu(Bc - diag(1./winv));
  X{s} = mldivide(U,mldivide(L,P*Bc*diag(winv./Sizes)));
  Xc_winv{s} = X{s}*Sizes - winv;
  ldCiBwi(s) = sum(log(winv).*(Sizes-1)) - sum(log(abs(diag(U))));
  oDinvo(s) = Sizes'*Xc_winv{s};  
end

kernel = struct('X',{X},...
		'oDinvo',{oDinvo},...
		'Xc_winv',{Xc_winv},...
	 	'ldCiBwi',{ldCiBwi});

if Data.Testing
  oDemes = Data.oDemes2;
  Jinvpt = Data.Jinvpt2;
  n = length(Jinvpt);
  o = ones(n,1);
  Jpt = sparse(1:n,Jinvpt,1);
  B = Ball(oDemes,oDemes);
  w = Wall(oDemes)*2;
  Delta = Jpt*B*Jpt' + Jpt*w*o'/2 + o*w'*Jpt'/2 - diag(Jpt*w);
  kernel.Delta = Delta;
end
