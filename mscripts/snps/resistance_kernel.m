

function kernel = resistance_kernel(Data,Graph,mValues,qValues)
%% Compute the auxiliary matrix X and the eigenvalues of X*inv(B) %%
%% necessary to compute efficiently the log likelihood terms      %%
%% Det(-inv(Dhat)*Q) and tr(-inv(Dhat)*Q*S)                       %%
%% This implementation uses the Schur complement trick described  %%
%% in (Hanks & Hooten, 2013) to avoid inverting a d-by-d matrix   %%
%% to compute all resistances when we only need the resistances   %%
%% between pairs of observed demes.                               %%
%% E. M. Hanks and M. B. Hooten. Circuit theory and model-based   %%
%% inference for landscape connectivity. Journal of the American  %%
%% Statistical	Association, 108(501):22-33, 2013                 %%
%% This implementation is more efficient for large graphs         %%


G = sparse(Graph.Vi,Graph.Vj,mValues);
L = diag(sum(G))-G; %% Graph Laplacian
J = Graph.Juniq;
K = Graph.Kuniq;

winv = reshape(1./qValues,[],1);

if Graph.allobsrv
  Hinv = L + 1;
else
  winv = winv(J);
  Hi11 = L(J,J) + 1;
  Hi22 = L(K,K) + 1;
  Hi21 = L(K,J) + 1;
  H22Hi21 = linsolve(Hi22,Hi21);
  Hinv = Hi11 - Hi21'*H22Hi21;
end
Binv = -2*Hinv;
if (Data.diploid)
  Binv = Binv/4;
  winv = winv/2;
end

Sizes = Data.Sizes;   %% the number of samples per deme
JtDJ = Data.JtDJ;     %% Jpt'*Dobs*Jpt
C = diag(Sizes); 

BinvW = Binv*diag(1./winv);
X = mldivide(C-BinvW,diag(winv));
ldCiBwi = logabsdet(X*BinvW);
cwinvt = Sizes*winv';
Bconst = winv'*JtDJ*winv;            %% winv'*J'*D*J*winv
oDinvoconst = - Sizes'*winv;         %% -diag(J*J')'*winv
ldDinvconst = log(sum(Sizes)) ...    %% log(n)
	    - sum(log(winv)) ...     %% logdet(diag(winv))
	    + sum(log(winv).*Sizes); %% logdet(diag(Jpt*winv))
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

if Data.Testing
  oDemes = Data.oDemes;
  Jinvpt = Data.Jinvpt;
  n = length(Jinvpt);
  o = ones(n,1);
  Jpt = sparse(1:n,Jinvpt,1);
  Ball = resistance_distance(G)/4;
  Wall = reshape(qValues,[],1);
  if (Data.diploid)
    Ball = 4*Ball;
    Wall = 2*Wall;
  end
  B = Ball(oDemes,oDemes);
  w = Wall(oDemes);
  %% Construct the entire n-ny-n fitted distance matrix
  Delta = Jpt*B*Jpt' + Jpt*w*o'/2 + o*w'*Jpt'/2 - diag(Jpt*w);
  kernel.Delta = Delta;
end
