

function X = AinvB(A,B,flag)


if (nargin==2)
  U = chol(A);
elseif ( ischar(flag) && strcmp(flag,'inv') )
  U = A;
end
Y = mldivide(U',B);
X = mldivide(U ,Y);
