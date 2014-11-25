

function ld = pseudologdet(A,r)


if nargin==1,r = rank(A);end
ld = sum(log(eigs(A,r)));
