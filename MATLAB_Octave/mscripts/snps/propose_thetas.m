

function [proposal,pi1_pi0] = propose_thetas(Data,kernel,params,schedule)

%%%%%%%%%%
type = 8;%
%%%%%%%%%%

thetai = schedule.paramtoupdate{type};
proposal.params = params;
proposal.type = type;
proposal.subtype = thetai;

n = Data.nIndiv;
df = params.df;

if (thetai==1)
  oDinvo = kernel.oDinvo;
  Xc_winv = kernel.Xc_winv;
  A = sum(sum(kernel.X.*Data.JtDJ));
  B = Xc_winv'*Data.JtDJ*Xc_winv;
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
