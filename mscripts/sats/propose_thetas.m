

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
    X = kernel.X{s};
    XC = kernel.XC{s};
    n = Sstruct.nIndiv(s);
    oDinvo = Sstruct.oDinvoconst(s) ...
           + sum(sum(X.*Sstruct.JtOJ{s}));
    A = sum(sum(X.*Sstruct.JtDJ{s}));
    B = Sstruct.Bconst(s) ...
      - sum(sum(X.*Sstruct.JtDJvct{s})) ...
      + sum(sum(XC'*Sstruct.JtDJ{s}*XC));
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
