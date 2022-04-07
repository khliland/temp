function [B, T, W, Q, R] = ncplsmem(ncomp, X, Y, Yadd)
% -----------------------------------------------------
% ----------- KH Liland + Ulf-mekk 2022 ---------------
% -----------------------------------------------------
% -- Solution of the Multiway Canonical PLS-problem ---
% -----------------------------------------------------
if nargin < 4
    Yadd = [];
end

% Initialize
nx = size(X); ndim  = length(nx);
ny = size(Y); nresp = ny(2);
W  = cell(1,ndim-1); wloads = W;
for i=1:(ndim-1)
    W{i} = zeros(nx(i+1), ncomp);
end
T = zeros(nx(1), ncomp); Pt  = double.empty([0,nx(2:end)]); % P = zeros([ncomp,nx(2:end)]);
Q = zeros(nresp, ncomp); Wa = double.empty([nx(2:end),0]); %Wa = zeros([nx(2:end),ncomp]); % 
ndWa = ndims(Wa);
% Centre the matrices:
X    = X - mean(X,1);       Y    = Y - mean(Y,1);
Yadd = Yadd - mean(Yadd,1); Yall = [Y, Yadd];

% Component loop
for a = 1:ncomp
    % (Candidate) loading weights
    W0 = GMPmem(GTmem(X,1), Yall, ndim-1); % W0 = X'*Yall
    if size(Yall,2) > 1               % Multiresponse and/or additional responses
        Z0 = GMPmem(X, W0,1);         % Z0 = X*W0
        [~,A] = ccXY(Z0,Y,[]);        % canoncorr(Z0,Y)
        w = GMPmem(W0,A(:,1),ndim-1); % w = W0*A
    else
        w = W0;                       % Single response, noe additional responses
    end
    
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
        end
    end
    
    % Apply sign convention
    for i = 1:length(wloads)
        sq = (wloads{i}.^2).*sign(wloads{i});
        wloads{i} = wloads{i}.*sign(sum(sq));
    end
    
    % Outerproducts of w
    wa = GOuter(wloads);
    
    % Scores
    t = GMPmem(X,wa,1); T(:,a) = t;
    
    % Loadings
    p = GMPmem(GTmem(X, 1), t, ndim-1)./(t'*t); % A matrix when X is 3d.
    q = Y'*t./(t'*t);  Q(:,a) = q;
    
    % Accumulate
    for i = 1:(ndim-1)
        W{i}(:,a) = wloads{i};
    end
    
    %P  = insert(P,p,a); %Wa = insert(Wa,wa,a,false); 
    Pt = cat(1,Pt,reshape(p, [1, size(p)]));  
    Wa = cat(ndWa,Wa,reshape(wa, [size(wa),1]));
    % Deflate
    X = X - GMPmem(t, reshape(p,[1,size(p)]), 1);
    Y = Y - t*q';
end
PtW = GMPmem(Pt,Wa,1); % Burde utnytte at PtW essensielt er øvre triangulær
R   = GMPmem(Wa, inv(PtW), ndim-1);
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

%% Outerproudct of vectors
function X = GOuter(x)
X = x{1};
if length(x) > 1
    for i=2:length(x)
        if i==2
            X = GMPmem(X, x{i}',1);
        else
            X = GMP(X, x{i},0); %OBS - hva skjer her i det generelle tilfellet?
        end
    end
end
