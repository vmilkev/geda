function A = scale( this, dat, col, flat, wndw )
% 
% scalind data to the range [0,1]
%
% INPUT:
%
% dat - data to scale, where
%
% col - column of dat() on which the scaling is based
%
% flat - adjust the data in a specific percentile to the same values (remove some peaks)
% wndw - is a specific percentile.
%
% OUTPUT:
%
% A(:,1) - time
% A(:,2) - value [0,1]

if ~exist('wndw','var') || isempty(wndw)
    wndw = 15;
end

right_cols =  dat(:,col) > 0.0 ;
dat2 = dat( right_cols,: ); % temporal solution; this should not be here

if ( isempty( dat2 ) )
    A = dat;
    return;
end

sz = size( dat2 );
A = zeros( sz(1,1),2 );

A(:,1) = dat2(:,1);

if ( flat ) % if we want to flattern data (cut some peacs)
    Ymax = prctile(dat2(:,col),100 - wndw);
    %Y50 = 0.5*max(dat2(:,col)) - 0.5*min(dat2(:,col));
    %Ymin = prctile(dat2(:,col),wndw);
    amax = dat2(:,col) > Ymax;
    amin = dat2(:,col) <= Ymax;
    maxmean = mean(dat2(amax,col));
    minmean = mean(dat2(amin,col));
    A(:,2) = dat2(:,col);

    a1 = A(amax,2)./maxmean - 1.0;
    a2 = A(amin,2)./minmean - 1.0;
    A(amax,2) = 1.0 + a1;
    A(amin,2) = -1.0 + a2;
    
else
%     Y = prctile(dat2(:,col),wndw);
%     A(:,2) = ( dat2(:,col) - min(dat2(:,col)) )./( max(dat2(:,col)) - min(dat2(:,col)) );
%     amin = A(:,2) <= Y;
%     A(amin,2) = A(amin,2) - 1.0;

    A(:,2) = dat2(:,col); % no scaling
end

