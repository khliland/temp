function [B, T, W, Q, R, Wa, Pt, X] = ncplsMissing(ncomp, X, Y, Yadd, orthW, outerW, impRep, impComp)
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
% impRep  - number imputation rounds
% impComp - number of components in imputations
if nargin < 4, Yadd = []; end
if nargin < 5, orthW = true; end
if nargin < 6, outerW = true; end
if nargin < 7, impRep = 10; end
if nargin < 8, impComp = ncomp; end

% Store for later
Xorig = X;
nvar  = ndims(X)-1;

% Initialize missing as 0
missing = ~isfinite(X);
X(missing) = normrnd(0,nanmean(X(:))/10,sum(missing(:)),1);

for i=1:impRep
    X = X-mean(X,1);
    [B, T, W, Q, R, Wa, Pt] = ncpls(ncomp, X, Y, Yadd, orthW, outerW);
    Xhat = GMP(GMP(X,R(:,:,1:impComp),nvar),Pt(1:impComp,:,:),1);
    X = Xorig;
    X(missing) = Xhat(missing);
end
