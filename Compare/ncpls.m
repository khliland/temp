function [B, T, W, Q, R, Wa, Pt] = ncpls(ncomp, X, Y, Yadd, orthW, outerW)
% -----------------------------------------------------
% ----------------- KH Liland 2022 --------------------
% -----------------------------------------------------
% -- Solution of the Multiway Canonical PLS-problem ---
% -----------------------------------------------------
% Input arguments
% ncomp   - number of components
% X       - data vector, matrix or tensor
% Y       - response vector or matrix
% Yadd    - (optional) addition response/meta information
% orthW   - orthogonalise loading weights on previous components (default = true)
% outerW  - use outerproduct of w_j (default = true) or full w
%           vector/matrix/tensor for deflation (unfolded version)
if nargin < 4, Yadd = []; end
if nargin < 5, orthW = true; end
if nargin < 6, outerW = true; end

% Initialize
nx = size(X); ndim  = length(nx);
ny = size(Y); nresp = ny(2);
W = cell(1,ndim-1); wloads = W;
for i=1:(ndim-1)
    W{i} = zeros(nx(i+1), ncomp);
end
T = zeros(nx(1), ncomp); Pt = zeros([ncomp,nx(2:end)]);
Q = zeros(nresp, ncomp); Wa = zeros([nx(2:end),ncomp]);

% Centre the matrices:
X    = X - mean(X,1);
Y    = Y - mean(Y,1);
Yadd = Yadd - mean(Yadd,1);

% Component loop
for a = 1:ncomp
    Yall = [Y, Yadd];
    
    % (Candidate) loading weights
    W0 = GMP(GT(X, ndim-1), Yall, 1); % W0 = X'*[Y Yadd]
    if size(Yall,2) > 1               % Multiresponse and/or additional response
        Z0 = GMP(X, W0, ndim-1);      % Z0 = X*W0
        if a > 1 && ~isempty(Yadd)    % Orthogonalise Z0 on previous components if Yadd is used
            Z0 = Z0 - T(:,1:(a-1))*(T(:,1:(a-1))'*Z0);
        end
        [~,A] = ccXY(Z0,Y,[]);        % canoncorr(Z0,Y)
        w = GMP(W0, A(:,1), 1);       % w = W0*A
    else
        w = W0;                       % Single response, noe additional response
    end
    
    % Choice of w vs outer(w_j) for deflation
    if outerW
        % Per-mode, normalized loading weights
        if ndim == 2 % Matrix X
            wloads{1} = w./norm(w);
        elseif ndim == 3 % Three-way X
            [w1,~,w2] = svd(w, 'econ');
            wloads{1} = w1(:,1);
            wloads{2} = w2(:,1);
        else % > three-way X
            wloads = parafac(w,1,[0 2 0 0 NaN]');
            for j = 1:length(wloads)
                wloads{j} = wloads{j}./norm(wloads{j});
                % Apply sign convention
                sq = (wloads{j}.^2).*sign(wloads{j});
                wloads{j} = wloads{j}.*sign(sum(sq)); 
            end
        end
    
        % Orthogonalise on previous components
        if orthW && a>1
            for i=1:(ndim-1)
                wloads{i} = wloads{i} - W{i}(:,1:(a-1))*(W{i}(:,1:(a-1))'*wloads{i});
                wloads{i} = wloads{i}./norm(wloads{i});
            end
        end

        % Outerproducts of w
        wa = GOuter(wloads);
    else
        % Orthogonalise on previous components
        if orthW && a>1
            Wvec = reshape(Wa,[prod(nx(2:end)),ncomp]);
            w = w(:) - Wvec*(Wvec'*w(:));
        end
        w = w(:)./norm(w(:));
        
        % Unfolded version
        if length(nx) > 2
            wa = reshape(w, nx(2:end));
        else
            wa = w;
        end

        % Dummies
        for i=1:(ndim-1)
            wloads{i} = zeros(nx(i+1),1);
        end
    end
    
    % Scores
    t = GMP(X,wa,ndim-1);
    if a>1
        t = t - T(:,1:(a-1))*(T(:,1:(a-1))'*t);
    end
    t = t./norm(t);
    
    % Loadings
    p = GMP(GT(X, ndim-1), t, 1); % Would be a matrix with 3D X
    q = Y'*t;
    
    % Accumulate
    T(:,a) = t;
    for i = 1:(ndim-1)
        W{i}(:,a) = wloads{i};
    end
    Q(:,a) = q;
    Pt  = insert(Pt,p,a);
    Wa = insert(Wa,wa,a,false);
    
    % Deflate
    Y = Y - t*q';
    if ~isempty(Yadd)
        Yadd = Yadd - t*(t'*Yadd);
    end
end
PtW = triu(GMP(Pt,Wa,ndim-1));
R   = GMP(Wa, eye(ncomp)/PtW, 1);
B   = cumsum(R.*reshape(Q',[ones(1,ndim-1),ncomp,nresp]),ndim);


%% Canonical correlations
function [r,A] = ccXY(X,Y,wt)
% Computes the coefficients in canonical variates between collumns of X and Y
[n,p1] = size(X); p2 = size(Y,2);

% Weighting of observations with regards to wt (asumes weighted centering already performed)
if ~isempty(wt)
    X = rscale(X,wt);
    Y = rscale(Y,wt);
end

% Factoring of data by QR decomposition and ellimination of internal linear
% dependencies in X and Y
[Q1,T11,perm1] = qr(X,0);       [Q2,T22,~] = qr(Y,0);
rankX          = sum(abs(diag(T11)) > eps(abs(T11(1)))*max(n,p1));
rankY          = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(n,p2));
if rankX < p1
    Q1 = Q1(:,1:rankX); T11 = T11(1:rankX,1:rankX);
end
if rankY < p2
    Q2 = Q2(:,1:rankY);
end

% Economical computation of canonical coefficients and canonical correlations
d = min(rankX,rankY);
if nargout == 1
    D    = svd(Q1' * Q2,0);
    r    = min(max(D(1:d), 0), 1); % Canonical correlations
else
    [L,D]    = svd(Q1' * Q2,0);
    A        = T11 \ L(:,1:d) * sqrt(n-1);
    % Transform back coefficients to full size and correct order
    A(perm1,:) = [A; zeros(p1-rankX,d)];
    r = min(max(diag(D(1:d)), 0), 1); % Canonical correlations
end


%% Insert into array
function Z = insert(Z,X,i,first)
ndim = length(size(Z));
if nargin < 4
    first = true;
end
if first
    switch ndim
        case 2
            Z(i,:) = X;
        case 3
            Z(i,:,:) = X;
        case 4
            Z(i,:,:,:) = X;
        case 5
            Z(i,:,:,:,:) = X;
        case 6
            Z(i,:,:,:,:,:) = X;
        case 7
            Z(i,:,:,:,:,:,:) = X;
        case 8
            Z(i,:,:,:,:,:,:,:) = X;
        case 9
            Z(i,:,:,:,:,:,:,:,:) = X;
        case 10
            Z(i,:,:,:,:,:,:,:,:,:) = X;
        otherwise
            error('Wrong size of Z')
    end
else % last
    switch ndim
        case 2
            Z(:,i) = X;
        case 3
            Z(:,:,i) = X;
        case 4
            Z(:,:,:,i) = X;
        case 5
            Z(:,:,:,:,i) = X;
        case 6
            Z(:,:,:,:,:,i) = X;
        case 7
            Z(:,:,:,:,:,:,i) = X;
        case 8
            Z(:,:,:,:,:,:,:,i) = X;
        case 9
            Z(:,:,:,:,:,:,:,:,i) = X;
        case 10
            Z(:,:,:,:,:,:,:,:,:,i) = X;
        otherwise
            error('Wrong size of Z')
    end
end


%% Outerproudct of vectors
function X = GOuter(x)
X = x{1};
if length(x) > 1
    for i=2:length(x)
        if i==2
            X = GMP(X, x{i}',1);
        else
            X = GMP(X, x{i},0);
        end
    end
end
