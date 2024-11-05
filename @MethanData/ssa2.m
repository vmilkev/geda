function [ H,l ] = ssa2( this, X, L, Q )

% X := signal to analyse
% L := window for Hankel matrix
% Q := number of eigen components to use in reconstruction, n <= L
%
% H := reconstructed signal
% l := vector of (all) eigenvalues, dim(l) = dim(L)

%--------------------------------------------------------------------

N = size(X,1);

K = N - L + 1;

t = (1:N)';

% trend = 0.001 * (t - 100).^2;
% p1 = 20;
% p2 = 30;
% periodic1 = 2 * sin(2*pi*t/p1);
% periodic2 = 0.75 * sin(2*pi*t/p2);
% 
% rng(123) % So we generate the same noisy time series every time.
% noise = 2 * (rand(N,1) - 0.5);
% X = trend(:,1) + periodic1(:,1) + periodic2(:,1) + noise(:,1);
% 
% figure(1);
% clf;
% hold on;
% plot(t,X);
% plot(t,trend);
% plot(t,periodic1);
% plot(t,periodic2);
% plot(t,noise);
% hold off;
%--------------------------------------------------------------------


% making Hankel matrix
H = zeros(L,K);
for i1 = 1:L
    H(i1,1:K) = X( i1 : i1 + K-1 );

end

% figure(2);
% imagesc(H);


% Hsqv = H*H';

% figure(3);
% imagesc(Hsqv);

[U,S,V] = svd(H);

% [coeff,score,latent,~,explained,mu] = pca(H);
% pca_95 = find(cumsum(explained)>50,1); % components explain more than 95% of all variability
% H95pca = coeff(:,1:pca_95)*coeff(:,1:pca_95)'*H;

noise_edge = 90; % procent of explained_eigenvalues

explained_svd = 100*diag(S.^2)./sum(diag(S.^2));
svd_expl = find(cumsum(explained_svd)>noise_edge,1); % components explain more than 95% of all variability

% Reconstructing
RC=zeros(N,L);
for m=1:L
  buf = S(m,m)*U(:,m)*V(:,m)';
  buf=buf(end:-1:1,:); % revert rows, in order to work with subdiagonals
  for n=1:N % anti-diagonal averaging
    ind = -(L-1)+n-1;
    RC(n,m)=mean( diag(buf,ind) ); % Get the elements on the (N-L+1)+n subdiagonal (k=-(N-L+1)+n) of A.
  end
end

% figure(5);
% subplot(1,2,1);
% plot(explained_svd);
% subplot(1,2,2);
% plot(cumsum(explained_svd));
% 
% figure(6);
% plot(t,X);
% hold on;
% for i = 1:12
%     plot(t,RC(:,i));
% end
% hold off;
% 
% figure(7);
% plot(t,X);
% hold on;
% plot(t,sum(RC(:,1:4),2));
% hold off;

w = zeros(N,1);
for i = 1:N
    if (i >= 1 && i <= L)
        w(i,1) = i + 1;
    elseif (i >= L+1 && i <= K)
            w(i,1) = L;
    elseif (i >= K+1 && i <= N)
            w(i,1) = N - i;
    end
end

W = zeros(L);
for i = 1:L
    for j = 1:L
        W(i,j) = ( RC(:,i)'*(RC(:,j).*w) )/( sqrt((RC(:,i)'*(RC(:,i).*w))) * sqrt((RC(:,j)'*(RC(:,j).*w))) );
    end
end

% figure(8);
% imagesc(W);

H = sum(RC(:,1:Q),2);
l = diag(S.^2);

end