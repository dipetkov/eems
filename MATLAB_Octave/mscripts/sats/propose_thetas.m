

function [proposal,pi1_pi0] = propose_thetas(Data,kernel,params,schedule)

%%%%%%%%%%
type = 8;%
%%%%%%%%%%

thetai = schedule.paramtoupdate{type};
proposal.params = params;
proposal.type = type;
proposal.subtype = thetai;

if (thetai==1)
  s2loc = params.s2loc;
  for s = 1:Data.nSites
    n = Data.nIndiv(s);
    oDinvo = kernel.oDinvo(s);
    Xc_winv = kernel.Xc_winv{s};
    A = sum(sum(kernel.X{s}.*Data.JtDJ{s}));
    B = Xc_winv'*Data.JtDJ{s}*Xc_winv;
    trDinvQxD = A - B/oDinvo;
    c = params.s2locShape + (n-1);
    d = params.s2locScale + trDinvQxD;
    s2loc(s) = rinvgam(c/2,d/2);
  end
  proposal.params.s2loc = s2loc;
  % This guarantees that the proposal is always accepted
  % i.e., a Gibbs update
  pi1_pi0 = Inf;
end
