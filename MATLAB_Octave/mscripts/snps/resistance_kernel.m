

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
J = Graph.Juniq;
K = Graph.Kuniq;
Hinv = diag(sum(G)) - G + 1;
winv = reshape(1./qValues,[],1);
JtDJ = Data.JtDJ;    %% Jpt'*Dobs*Jpt
Sizes = Data.Sizes;  %% the number of samples per deme

if Graph.allobsrv
  Binv = - 0.5 * Hinv;
else
  winv = winv(J);
  Binv = - 0.5 * Hinv(J,J) + 0.5 * Hinv(J,K) * linsolve( Hinv(K,K) , Hinv(K,J) );
end
if (Data.diploid)
  winv = winv/2;
else
  Binv = Binv*4;
end

Bc = linsolve(Binv,diag(Sizes));
[L,U,P] = lu(Bc - diag(1./winv));
X = mldivide(U,mldivide(L,P*Bc*diag(winv./Sizes)));
Xc_winv = X*Sizes - winv;
ldCiBwi = sum(log(winv).*(Sizes-1)) - sum(log(abs(diag(U))));
oDinvo = Sizes'*Xc_winv;

kernel = struct('X',{X},...
		'oDinvo',{oDinvo},...
		'Xc_winv',{Xc_winv},...
	 	'ldCiBwi',{ldCiBwi});

if Data.Testing
  oDemes = Data.oDemes;
  Jinvpt = Data.Jinvpt;
  n = length(Jinvpt);
  o = ones(n,1);
  Jpt = sparse(1:n,Jinvpt,1);
  Ball = resistance_distance(G);
  Wall = reshape(qValues,[],1);
  if (Data.diploid)
    Wall = Wall*2;
  else
    Ball = Ball/4;
  end
  B = Ball(oDemes,oDemes);
  w = Wall(oDemes);
  %% Construct the entire n-ny-n fitted distance matrix
  Delta = Jpt*B*Jpt' + Jpt*w*o'/2 + o*w'*Jpt'/2 - diag(Jpt*w);
  kernel.Delta = Delta;
end
