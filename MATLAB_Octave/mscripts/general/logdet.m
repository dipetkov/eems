function y = logdet(A,flag)
% log(det(A)) where A is positive-definite.
% This is faster and more stable than using log(det(A)).

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.


if (nargin==1)
  U = chol(A);
elseif ( ischar(flag) && strcmp(flag,'inv') )
  U = A;
end

y = 2*sum(log(diag(U)));
