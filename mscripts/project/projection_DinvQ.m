

function DinvQ = projection_DinvQ(D)


n = nrow(D);
onev = ones(n,1);
Dinvo = mldivide(D,onev);
oDinvo = sum(Dinvo);
Q = eye(n) - onev*Dinvo'/oDinvo;
DinvQ = mldivide(D,Q);
