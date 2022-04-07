%% Generalised orthogonalisation, 
% Orthgonoalise one tensor on another
function X = Gorth(X, Y, Yorth)
% X       - tensor to be orthogonalised
% Y       - fixed tensor
% Yorth   - if true, assume Y is orthonormal
% nactive - number of right facing dimensions, e.g. active dimensions.


% Dimensions
nx = size(X);
nactive = length(nx);
ny = size(Y);
if length(ny) == nactive
    ny = [ny,1];
end

X = reshape(X,[prod(nx),1]);
Y = reshape(Y,[prod(ny(1:nactive)), prod(ny((nactive+1):end))]);
if Yorth
    X = reshape(X - Y*(Y'*X), nx);
else
    X = reshape(X - (Y/(Y'*Y))*(Y'*X), nx);
end