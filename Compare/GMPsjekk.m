% {
clear, clc
n = 10; m = 5; p = 7; q = 3; s = 8; v = 4; t = 4;
A = randn(n,1); B = randn(n,m);
C = randn(n,p); D = randn(n,q,s);
E = randn(q,s); F = randn(q,v,t);
G = randn()
vec = @(a)a(:); dev = inf(11,1);

% (1xn) x (nx1)
AAm = GMPmem(A', A,1);
AA = GMP(A', A);
dev(1) = max(abs(vec(AA - A'*A))); % Check acuracy

% (nx1) x (1xn)
AAtm = GMPmem(A, A',1);
AAt = GMP(A, A');
dev(2) = max(abs(vec(AAt - A*A'))); % Check acuracy

% (mxn) x (nxm)
BBm = GMPmem(B', B,1);
BB = GMP(B', B);
dev(3) = max(abs(vec(BB - B'*B))); % Check acuracy

% (nxm) x (mxn)
BBtm = GMPmem(B, B',1);
BBt = GMP(B, B');
dev(4) = max(abs(vec(BBt - B*B'))); % Check acuracy

% (1xn) x (nxm)
ABm = GMPmem(A', B,1);
AB = GMP(A', B);
dev(5) = max(abs(vec(AB - A'*B))); % Check acuracy

% (mxn) x (nxp)
BCm = GMPmem(B', C,1);
BC = GMP(B', C);
dev(6) = max(abs(vec(BC - B'*C))); % Check acuracy
norm(BC-B'*C, 'fro')

% (1xn) x (nxqxs)
ADm = GMPmem(A', D,1);
AD = GMP(A', D);
dev(7) = sum([1,q,s] ~= size(AD)); % Check dimensions

% (nxqxs) x (qxs) i.e. (nxqxs) x (qxsx1) 
DEm = GMPmem(D, E, 1);
DE = GMP(D, E, 2);
dev(8) = sum([n,1] ~= size(DE)); % Check dimensions

% ------------------------------
% (nxqxs) x (qxsxn) : qxs inside
[Dt, dactive] = GT(D, 2); % Generalised transpose
DD1m = GMPmem(D, Dt, 1);
DD1 = GMP(D, Dt, 2);
dev(9) = sum([n,n] ~= size(DD1)); % Check dimensions

% (qxsxn) x (nxqxs) : n inside
DD2m = GMPmem(Dt, D, 2);
DD2 = GMP(Dt, D, dactive);
dev(10) = sum([q,s,q,s] ~= size(DD2)); % Check dimensions

% ------------------------------
% Three arrays : when all nactive == 1 it seems like (A*B)*C == A*(B*C)
Et = GTmem(E,1);
DEF1m = GMPmem(GMPmem(D,Et,2),F,2);
DEF2m = GMPmem(D,GMPmem(Et,F,1),2);
dev(11) = max(abs(DEF1m(:)-DEF2m(:)));

Et = GT(E);
DEF1 = GMP(GMP(D,Et),F);
DEF2 = GMP(D,GMP(Et,F));
dev(11) = max(abs(DEF1(:)-DEF2(:)));

disp(dev)
% Three arrays : when nactive > 1, (A*B)*C ~= A*(B*C)
% DEF3m = GMP(GMP(D,E,1),F,2);
% DEF3 = GMP(GMP(D,E,2),F);
% DEF4m = GMP(D,GMP(E,F,2),1);
% DEF4 = GMP(D,GMP(E,F),2); % <- Wrong dimension since second dim of E is lost
%}
