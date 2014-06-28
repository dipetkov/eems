

function [proposal,pi1_pi0] = propose_thetas(Sstruct,kernel,params,schedule)

%%%%%%%%%%
type = 1;%
%%%%%%%%%%

thetai = schedule.paramtoupdate{type};
proposal.params = params;
proposal.type = type;
proposal.subtype = thetai;

if (thetai==1)
  s2loc = params.s2loc;
  for s = 1:Sstruct.nSites
    n = Sstruct.nIndiv(s);
    trDinvQxD = params.trDinvQxD(s);
    c = params.s2locC + (n-1);
    d = params.s2locD + trDinvQxD;
    s2loc(s) = rinvgam(c/2,d/2);
  end
  proposal.params.s2loc = s2loc;
  % This guarantees that the proposal is always accepted
  % i.e., a Gibbs update
  pi1_pi0 = Inf;
end
