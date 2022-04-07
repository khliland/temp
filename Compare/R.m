    function [w, s] = R(X, W, Yp, t)
        % XW - centred and transformed X-data matrix(independant), Yp - matrix of dependent variables, W - temporary loading weights
        XW = X*W; if nargin > 3, XW = XW - t*(t'*XW); end 
        [A,~,r] = ccXY(XW, Yp); % Computaion of canonical correlations between XW and Y.
        w       = W*A(:,1); w = w/norm(w); % The optimal loading weight vector of norm 1
        s       = r(1)^2; % squared canonical correlation