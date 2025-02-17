function res = apply_ssa_ams(this, arr, max_length, ssa_window, ssa_eigenvalue)

res = zeros(size(arr));
Acell = {};
n_waves = floor( numel(arr)/max_length ); % number of waves (sub-signals)
r1 = zeros(n_waves);
r2 = zeros(n_waves);
r1(1) = 1;
r2(1) = max_length;
for ii = 2:n_waves
    r1(ii) = r1(ii-1) + max_length;
    r2(ii) = r2(ii-1) + max_length;
end
r2(n_waves) = numel(arr);

parfor ii = 1:n_waves
    %[A0,~] = this.ssa2(arr(r1(ii):r2(ii)),ssa_window,ssa_eigenvalue);
    [A0,~,~,~] = in_ssa(arr(r1(ii):r2(ii)),ssa_window,ssa_eigenvalue,[]); % this way is 4 times faster!
    Acell{ii} = A0;
end
for ii = 1:n_waves
    res( r1(ii):r2(ii) ) = Acell{ii};
end

end

%%

function [ H,l,nl,F ] = in_ssa( X, L, Q, E )

% X := signal to analyse
% L := size of window for Hankel matrix
% Q := number of eigen components to use in reconstruction, n <= L;
%      if Q == [], Q will be defined based on argument E
% E := procent of variance explained by eigenvalues;
%      this value defines number of eigenvalues (Q)
%      used for reconstruction; this option is active if Q == []
%
% H := reconstructed signal
% l := vector of (all) eigenvalues, dim(l) = dim(L)
% nl := number of eigenvalues responsible for E % variation
% F := first nl number of eigenvectors

N = size(X,1);
t = (1:N)';

% Calculate covariance matrix C (Toeplitz approach)

% Calculate the covariance matrix.
% There are several numerical approaches to estimate C.
% Here, we calculate the covariance function with CORR and build C with the function TOEPLITZ.

% covX = xcorr(X,L-1,'unbiased');
% Ctoep=toeplitz(covX(L:end));

% Calculate covariance matrix (trajectory approach)

% An alternative approach is to determine C directly from the scalar product of Y,
% the time-delayed embedding of X.
% Although this estimation of C does not give a Toeplitz structure,
% with the eigenvectors not being symmetric or antisymmetric,
% it ensures a positive semi-definite covariance matrix.

Y=zeros(N-L+1,L);
for m=1:L
    Y(:,m) = X((1:N-L+1)+m-1);
end

Cemb=Y'*Y / (N-L+1);

% Choose covariance estimation
% Choose between Toeplitz approach (cf. Vautard & Ghil)
% and trajectory approach (cf. Broomhead & King).

% C=Ctoep;
C=Cemb;

% Calculate eigenvalues LAMBDA and eigenvectors V

% In order to determine the eigenvalues and eigenvectors of C,
% we use the function EIG. This function returns two matrices,
% the matrix V with eigenvectors arranged in columns,
% and the matrix LAMBDA with eigenvalues along the diagonal.

[V,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);               % extract the diagonal elements
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
V = V(:,ind);                    % and eigenvectors

% Calculate principal components PC
% The principal components are given as the scalar product between Y,
% the time-delayed embedding of X, and the eigenvectors V.

PC = Y*V;

% Calculate reconstructed components RC
% In order to determine the reconstructed components RC,
% we have to invert the projecting PC = Y*V;
% i.e. RC = Y*V*V'=PC*V'.
% NOTE, V is orthogonal matrix, therefore inv(V) = V'.
% Averaging along anti-diagonals gives the RCs for the original input X.

RC=zeros(N,L);
for m=1:L
    buf=PC(:,m)*V(:,m)'; % invert projection
    buf=buf(end:-1:1,:); % revert rows, in order to work with subdiagonals
    for n=1:N % anti-diagonal averaging
        RC(n,m)=mean( diag(buf,-(N-L+1)+n) ); % Get the elements on the (N-L+1)+n subdiagonal (k=-(N-L+1)+n) of A.
    end
end

% Compare reconstruction and original time series
if ~exist('E','var') || isempty(E)
    if isempty(Q)
        disp('Error: At least one of these two parameters should be supplied: Q or E!');
        return
    else
        H = sum(RC(:,1:Q),2);
        nl = Q;
    end
else
    l_sum = cumsum(LAMBDA)*100./sum(LAMBDA);
    n_E = find( l_sum <= E );
    if isempty(n_E)
        nl = 1;
    else
        nl = max(n_E) + 1;
    end
    H = sum(RC(:,1:nl),2);
end

l = LAMBDA;

F = V;%(:,1:nl);

end
