

function kernel = resistance_kernel(Sstruct,Mstruct,mValues)
%% Compute the auxiliary matrix X and the eigenvalues of X*inv(B) %%
%% necessary to compute efficiently the log likelihood terms      %%
%% Det(-inv(Delta)*Q) and tr(-inv(Delta)*Q*S)                     %%


G = sparse(Mstruct.Mi,Mstruct.Mj,mValues);
Ball = resistance_distance(G)/4;
Wall = ones(Mstruct.nDemes,1);

if (Sstruct.diploid)
   Ball = 4*Ball;
   Wall = 2*Wall;
end

Sizes = Sstruct.Sizes;   %% the number of samples per deme
oDemes = Sstruct.oDemes; %% the observed demes
Cinv = diag(1./Sizes);
C = diag(Sizes); 
B = Ball(oDemes,oDemes);
w = Wall(oDemes);
BWinv = B*diag(1./w);
X = mldivide(B*C-diag(w),BWinv);
ldCiBwi = logabsdet(Cinv-BWinv);
kernel = struct('X',{X},'XC',{X*C},...
	 	'ldCiBwi',{ldCiBwi});
