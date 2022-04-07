function [B, W, P, T, mX, cc] = cpls(pc, X, Yprim, Ysec)
% Function for estimating PLS-loading weights based on canonical correlation between the 
% transformed data XV, where V = X'*Y with Y = [Yprim Ysec] contain the directions 
% maximizing the indivudal covariances between X and each column of Y.
% ---------------
% INPUTS:
% ---------------
% pc    - number of components to be extracted
% X     - matrix of independent variables
% Yprim - matrix of primary dependent variables
% Ysec  - matrix of secondary dependent variables
% wt    - weighting of the observations (rows in X and Y)
% ---------------
% OUTPUTS:
% ---------------
% B  - regression coeffs
% Q  - loading weights
% P  - P-loadings
% T  - scores
% mX - column means of X
% cc - squared canonical correlations for the components
% ---------------
[n,m] = size(X); nyp = size(Yprim,2); %nys = size(Ysec,2);
mX = mean(X); X = X - mX; % Centering in case of no weighting.
if nargin < 4, Y = Yprim; else Y = [Yprim Ysec]; end
mY = mean(Y); Y = Y - mY;
W  = zeros(m,pc); P = zeros(m,pc);
T  = zeros(n,pc); Q  = zeros(pc,nyp);
cc = zeros(pc,1);

for j=1:pc  
   W0 = X'*Y; %[Yprim Ysec]; 
   if j == 1, [w, cc(j)] = R(X, W0, Y(:,1:nyp)); % Calculation of loading weights by canonical correlation
   else [w, cc(j)] = R(X, W0, Y(:,1:nyp),T(:,1:j-1)); end %, w = w - W(:,1:j-1)*(W(:,1:j-1)'*w); w = w/norm(w);
   t = X*w; if j > 1, t = t - T(:,1:j-1)*(T(:,1:j-1)'*t); end
   t = t/norm(t); p = X'*t; q = t'*Y;
   Y = Y - t*q;  
   W(:,j) = w; T(:,j) = t; P(:,j) = p; Q(j,:) = q(1:nyp);
end
% B = cumsum(W/triu(P'*W).*Q(:,1)',2);
% B = [mY(1:nyp)-mX*B; B];
B = regcoeffs(W, triu(P'*W), Q, mX, mY(1:nyp)); % the regressions coeffs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Internal functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [w, s] = R(X, W, Yp, t)
        % XW - centred and transformed X-data matrix(independant), Yp - matrix of dependent variables, W - temporary loading weights
        XW = X*W; if nargin > 3, XW = XW - t*(t'*XW); end 
        [A,~,r] = ccXY(XW, Yp); % Computaion of canonical correlations between XW and Y.
        w       = W*A(:,1); w = w/norm(w); % The optimal loading weight vector of norm 1
        s       = r(1)^2; % squared canonical correlation
%------------

    function [A,B,r] = ccXY(X,Y) %,wt)
        % Function that finds the loading weights (ceoefficients) of the canonical variates of X and Y
        % This is a "stripped" version of the function canoncorr.m from the Statistics TB.
        [n,p1] = size(X); p2 = size(Y,2);
        % Factoring of data by QR-decomposition and elimination of internal linear
        % dependencies in X and Y, respectively:
        [Q1,T11,perm1] = qr(X,0);       [Q2,T22,perm2] = qr(Y,0);
        rankX          = sum(abs(diag(T11)) > eps(abs(T11(1)))*max(n,p1));
        rankY          = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(n,p2));
        if rankX < p1
            %     warning('canoncorr:NotFullRank','X is not full rank.');
            Q1 = Q1(:,1:rankX); T11 = T11(1:rankX,1:rankX);
        end
        if rankY < p2
            %     warning('canoncorr:NotFullRank','Y is not full rank.');
            Q2 = Q2(:,1:rankY); T22 = T22(1:rankY,1:rankY);
        end
        % Economical estimation of the canonical coefficients and the canonical correlations:
        d = min(rankX,rankY);
        [L,D,M]    = svd(Q1' * Q2,0);
        A          = T11 \ L(:,1:d) * sqrt(n-1);
        B          = T22 \ M(:,1:d) * sqrt(n-1);
        r          = min(max(diag(D(:,1:d))', 0), 1); % Canonical correlations with elimination of roundoff errors
        % Backtransformation of coefficients to full size and correct order:
        A(perm1,:) = [A; zeros(p1-rankX,d)];
        B(perm2,:) = [B; zeros(p2-rankY,d)];
%------------

    function betas = regcoeffs(W, R, qy, mX, mY)
        [nvar,q] = size(qy);
        betas = cumsum(repmat(W/R,1,1,q).*reshape(qy, 1, nvar,q),2); % eye(nvar)/R; - the inverse of R by "backsubstitution".
        betas = cat(1, reshape(mY,1,1,[])-sum(mX'.*betas), betas); % Append the constant terms for all models
            
