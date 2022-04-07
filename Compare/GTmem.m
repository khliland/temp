% Generalised transpose 
function Xt = GTmem(X, ox)
% INPUT:
% X - array 
% ox - the number of outer x-dimensions shifted in the transposition.
% OUTPUT:
% Xt - the generalized transpose
perm = [ox+1:ndims(X) 1:ox];
Xt   = permute(X,perm);