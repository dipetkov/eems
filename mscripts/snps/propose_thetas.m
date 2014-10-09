

function [proposal,pi1_pi0] = propose_thetas(Data,kernel,params,schedule)

%%%%%%%%%%
type = 1;%
%%%%%%%%%%

thetai = schedule.paramtoupdate{type};
proposal.params = params;
proposal.type = type;
proposal.subtype = thetai;

n = Data.nIndiv;
df = params.df;

if (thetai==1)
  XC = kernel.XC;
  X = kernel.X;
  oDinvo = kernel.oDinvoconst ...
         + sum(sum(X.*Data.JtOJ));
  A = sum(sum(X.*Data.JtDJ));
  B = kernel.Bconst ...
    - sum(sum(X.*kernel.cvtJtDJvct)) ...
    + sum(sum(XC'*Data.JtDJ*XC));
  trDinvQxD = A - B/oDinvo;
  c = params.s2locShape + df*(n-1);
  d = params.s2locScale + df*trDinvQxD;
  s2loc = rinvgam(c/2,d/2);
  proposal.params.s2loc = s2loc;
  % This guarantees that the proposal is always accepted
  % i.e., a Gibbs update
  pi1_pi0 = Inf;
elseif (thetai==2)
  df = rnorm(params.df,params.dfProposalS2);
  if (df>params.dfmin && df<params.dfmax)
    % This corresponds to pi(df)\propto(1/df)
    pi1_pi0 = log(params.df)-log(df);
  else
    pi1_pi0 = -Inf;
  end
  proposal.params.df = df;
end
