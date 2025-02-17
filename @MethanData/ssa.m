function [ H,l,nl,F ] = ssa( this, X, L, Q, E )

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

%% Calculate covariance matrix C (Toeplitz approach)

% Calculate the covariance matrix.
% There are several numerical approaches to estimate C.
% Here, we calculate the covariance function with CORR and build C with the function TOEPLITZ.

% covX = xcorr(X,L-1,'unbiased');
% Ctoep=toeplitz(covX(L:end));

% figure(2);
% set(gcf,'name','Covariance matrix');
% clf;
% imagesc(Ctoep);
% axis square
% set(gca,'clim',[-1 1]);
% colorbar

%% Calculate covariance matrix (trajectory approach)

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

% figure(22);
% set(gcf,'name','Covariance matrix');
% clf;
% imagesc(Cemb);
% axis square
% set(gca,'clim',[-1 1]);
% colorbar

%% Choose covariance estimation

% Choose between Toeplitz approach (cf. Vautard & Ghil)
% and trajectory approach (cf. Broomhead & King).

% C=Ctoep;
C=Cemb;

%% Calculate eigenvalues LAMBDA and eigenvectors V

% In order to determine the eigenvalues and eigenvectors of C,
% we use the function EIG. This function returns two matrices,
% the matrix V with eigenvectors arranged in columns,
% and the matrix LAMBDA with eigenvalues along the diagonal.

[V,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);               % extract the diagonal elements
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
V = V(:,ind);                    % and eigenvectors

% figure(3);
% set(gcf,'name','Eigenvectors V and eigenvalues LAMBDA')
% clf;
% subplot(3,1,1);
% plot(LAMBDA,'o-');
% subplot(3,1,2);
% plot(V(:,1:2), '-');
% legend('1', '2');
% subplot(3,1,3);
% plot(V(:,3:4), '-');
% legend('3', '4');


%% Calculate principal components PC

% The principal components are given as the scalar product between Y,
% the time-delayed embedding of X, and the eigenvectors V.

PC = Y*V;
%Y2 = PC*V';

% figure(4);
% set(gcf,'name','Principal components PCs')
% clf;
% for m=1:4
%   subplot(4,1,m);
%   plot(t(1:N-L+1),PC(:,m),'k-');
%   ylabel(sprintf('PC %d',m));
%   ylim([-10 10]);
% end

%% Calculate reconstructed components RC

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

% figure(5);
% set(gcf,'name','Reconstructed components RCs')
% clf;
% components = 7;
% for m=1:components
%   subplot(components,1,m);
%   plot(t,RC(:,m),'r-');
%   ylabel(sprintf('RC %d',m));
%   ylim([-1 1]);
% end

%% Compare reconstruction and original time series

% Note that the original time series X can be completely reconstructed
% by the sum of all reconstructed components RC (upper panel).
% The sine function can be reconstructed with the first pair of RCs (lower panel).

% figure(6);
% set(gcf,'name','Original time series X and reconstruction RC')
% clf;
% subplot(2,1,1)
% plot(t,X,'b.',t,sum(RC(:,:),2),'r-');
% legend('Original','Complete reconstruction');
% 
% subplot(2,1,2)
% %plot(t,X,'b.','LineWidth',1);
% plot(t,X,'b.-.',t,sum(RC(:,1:2),2),'r-');
% %plot(t, X0, 'b.-',t,sum(RC(:,1:2),2),'r-');
% %legend('Original','Reconstruction with RCs 1-2');

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