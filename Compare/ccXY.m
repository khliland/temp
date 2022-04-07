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