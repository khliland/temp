function Z = GMPmem(X, Y, ox)
% Generalized matrix product function
% INPUT: 
% X and Y - arrays of matching dimensions .
% ox - the number of outer x-dimensions in the multiplication.
% OUTPUT: 
% Z - the generalized matrix product of X and Y.
sy = [size(Y) 1]; sx = size(X); 
ix = ndims(X)-ox;                          % ix - the number of inner x(&y)-dimensions in the multiplication.
X  = reshape(X, prod(sx(1:ox))    ,[]);    % Convert X-array to matrix
Y  = reshape(Y, prod(sx(ox+1:end)),[]);    % Convert Y-array to matrix
Z  = reshape(X*Y,[sx(1:ox),sy(ix+1:end)]); % Multiply and reshape to correct dims.
end
