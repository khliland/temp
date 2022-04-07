%% Generalised transpose, 
% Switch active and passive dimensions of an array
function [X, nleft] = GT(X, nright)
% nright - number of right facing dimensions, e.g. active dimensions.
%        - defaults to 1 if not specified
if nargin == 1
    nright = 1;
end
nleft = length(size(X))-nright;
perm  = [(nleft+1):length(size(X)) (1:nleft)];
X = permute(X, perm);
