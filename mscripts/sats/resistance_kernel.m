

function kernel = resistance_kernel(Sstruct,Mstruct,mValues)
%% Compute the auxiliary matrix X and the eigenvalues of X*inv(B) %%
%% necessary to compute efficiently the log likelihood terms      %%
%% Det(-inv(Delta)*Q) and tr(-inv(Delta)*Q*S)                     %%


G = sparse(Mstruct.Mi,Mstruct.Mj,mValues);
Ball = resistance_distance(G)/4;
Wall = ones(Mstruct.nDemes,1);

%if (Sstruct.microsat)
%   Ball = Ball/4;
%   Wall = Wall/2;
%end

nSites = Sstruct.nSites;
X = cell(nSites,1);
XC = cell(nSites,1);
ldCiBwi = zeros(nSites,1);

for s = 1:nSites
  Sizes = Sstruct.Sizes{s};   %% the number of samples per deme
  oDemes = Sstruct.oDemes{s}; %% the observed demes
  Cinv = diag(1./Sizes);
  C = diag(Sizes); 
  B = Ball(oDemes,oDemes);
  w = Wall(oDemes);
  BWinv = B*diag(1./w);
  X{s} = mldivide(B*C-diag(w),BWinv);
  XC{s} = X{s}*C;
  ldCiBwi(s) = logabsdet(Cinv-BWinv);
end

kernel = struct('X',{X},'XC',{XC},...
	 	'ldCiBwi',{ldCiBwi});
