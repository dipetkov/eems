

function kernel = resistance_kernel(Sstruct,Mstruct,mValues,qValues)
%% Compute the auxiliary matrix X and the eigenvalues of X*inv(B) %%
%% necessary to compute efficiently the log likelihood terms      %%
%% Det(-inv(Dhat)*Q) and tr(-inv(Dhat)*Q*S)                       %%


G = sparse(Mstruct.Mi,Mstruct.Mj,mValues);
Ball = resistance_distance(G)/4;
Wall = reshape(qValues,[],1);

if (Sstruct.diploid)
   Ball = 4*Ball;
   Wall = 2*Wall;
end

Sizes = Sstruct.Sizes;   %% the number of samples per deme
oDemes = Sstruct.oDemes; %% the observed demes
JtDJ = Sstruct.JtDJ;
Cinv = diag(1./Sizes);
C = diag(Sizes); 
B = Ball(oDemes,oDemes);
w = Wall(oDemes);
winv = 1./w;
BWinv = B*diag(winv);
X = mldivide(B*C-diag(w),BWinv);
ldCiBwi = logabsdet(Cinv-BWinv);
cwinvt = Sizes*winv';
Bconst = winv'*JtDJ*winv;            %% winv'*J'*D*J*winv
oDinvoconst = - Sizes'*winv;         %% -diag(J*J')'*winv
ldDinvconst = log(sum(Sizes)) ...
	    - sum(log(Sizes)) ...    %% logdet(C) = logdet(diag(Sizes))
	    + sum(log(winv).*Sizes); %% logdet(W) = logdet(diag(J*w))
cvtJtDJvct = JtDJ*cwinvt' + cwinvt*JtDJ';

kernel = struct('X',{X},'XC',{X*C},...
		'Bconst',{Bconst},...
	 	'ldCiBwi',{ldCiBwi},...
		'cvtJtDJvct',...
		{cvtJtDJvct},...
		'ldDinvconst',...
		{ldDinvconst},...
		'oDinvoconst',...
		{oDinvoconst});
