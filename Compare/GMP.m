%{
n = 10; m = 5; p = 7; q = 3; s = 8; v = 4;
A = randn(n,1);
B = randn(n,m);
C = randn(n,p);
D = randn(n,q,s);
E = randn(q,s);
F = randn(q,v);
G = randn(q,s,v);
vec = @(a)a(:);
dev = inf(1,10);

% (1xn) x (nx1)
AA = GMP(A', A);
dev(1) = max(abs(vec(AA - A'*A))); % Check acuracy

% (nx1) x (1xn)
AAt = GMP(A, A');
dev(2) = max(abs(vec(AAt - A*A'))); % Check acuracy

% (mxn) x (nxm)
BB = GMP(B', B);
dev(3) = max(abs(vec(BB - B'*B))); % Check acuracy

% (nxm) x (mxn)
BBt = GMP(B, B');
dev(4) = max(abs(vec(BBt - B*B'))); % Check acuracy

% (1xn) x (nxm)
AB = GMP(A', B);
dev(5) = max(abs(vec(AB - A'*B))); % Check acuracy

% (mxn) x (nxp)
BC = GMP(B', C);
dev(6) = max(abs(vec(BC - B'*C))); % Check acuracy
norm(BC-B'*C, 'fro')

% (1xn) x (nxqxs)
AD = GMP(A', D);
dev(7) = sum([1,q,s] ~= size(AD)); % Check dimensions

% (nxqxs) x (qxs) i.e. (nxqxs) x (qxsx1) 
DE = GMP(D, E, 2);
dev(8) = sum([n,1] ~= size(DE)); % Check dimensions

% (nxqxs) x (qxsxn) : qxs inside
[Dt, dactive] = GT(D, 2); % Generalised transpose
DD1 = GMP(D, Dt, 2);
dev(9) = sum([n,n] ~= size(DD1)); % Check dimensions

% (qxsxn) x (nxqxs) : n inside
DD2 = GMP(Dt, D, dactive);
dev(10) = sum([q,s,q,s] ~= size(DD2)); % Check dimensions

% Three arrays : when all nactive == 1 it seems like (A*B)*C == A*(B*C)
Et = GT(E);
DEF1 = GMP(GMP(D,Et),F);
DEF2 = GMP(D,GMP(Et,F));
dev(11) = max(abs(DEF1(:)-DEF2(:)));

% Three arrays : when nactive > 1, (A*B)*C ~= A*(B*C)
% DEF3 = GMP(GMP(D,E,2),G);
% DEF4 = GMP(D,GMP(E,G),2); % <- Wrong dimension since second dim of E is lost

disp(dev)
%}

%% Generalized matrix product, equivalent to https://numpy.org/doc/stable/reference/generated/numpy.tensordot.html
function Z = GMP(X, Y, nactive)
% nactive - number of active dimensions of X and Y (shared dimensions)
%         - defaults to 1 if not specified
%         - if equal to 0 it results in "element-wise outer-product" 
%         - when nactive > 1, (A*B)*C ~= A*(B*C)
%         - when all nactive == 1 it seems like (A*B)*C == A*(B*C)
if nargin < 3
    nactive = 1;
end

% Dimensions
nx = size(X);
ny = size(Y);
xout = length(nx) - nactive; % Number of passive dimensions for X
if length(ny) == nactive
    ny = [ny,1];
end

Z = reshape( ...
    reshape(X,[prod(nx(1:xout)),prod(ny(1:nactive))]) * ...
    reshape(Y,[prod(ny(1:nactive)), prod(ny((nactive+1):end))]), ...
    [nx(1:xout), ny((nactive+1):end)]);
end