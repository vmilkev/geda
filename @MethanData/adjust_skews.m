function res = adjust_skews( this, adj_table, lim_range, lim_time )

% adj_table - alignment results matrix
% lim_range - number of records in AMS signal
% lim_time - emission signal length, sec 

s = adj_table(:,1); % estimated time skews

d = fnd_clst( s, lim_time );

if ( ~any(d{1}) )
    res = [];
    disp('There are no clusters found in the data to be adjusted.');
    return;
end

res( 1:size(s,1), 1:8 ) = 0;

for i = 1:size(d,2)
    j_correct = d{i}(:,2);
    res( j_correct, 1 ) = d{i}(:,1); % correct timing skews
    res( j_correct, 2 ) = 1; % flags for correct timing skews
    % ---------------------------
    res(j_correct,3) = adj_table(j_correct,1);
    res(j_correct,4) = adj_table(j_correct,2);
    res(j_correct,5) = adj_table(j_correct,3);
    % ---------------------------
    res(:,6) = adj_table(:,4);
    res(:,7) = adj_table(:,5);
    res(:,8) = adj_table(:,6);
    % ---------------------------
end

i = 1;
j = i;

while i <= size(s,1)
    if ( res(i,2) )
        range = i - j;
        if ( ~range )
            j = i + 1;
            i = i + 1;
            continue;
        else
            n = i-1;
            while n >= j
                res(n,1) = res(i,1);
                res(n,3) = adj_table(i,1);
                res(n,5) = res(n+1,4);
                res(n,4) = res(n,5) - ( res(n,7) - res(n,6) );
                if ( res(n,4) > lim_range )
                    res(n,4) = 0;
                end
                if ( res(n,5) > lim_range )
                    res(n,5) = 0;
                end
                n = n - 1;
            end
            j = i+1;
        end
    end
    i = i + 1;
    if ( i == size(s,1) && ~res(i,2) )
        ii = i;
        while ~res(ii,2)
            ii = ii - 1;
        end
        n = ii+1;
        while n <= i
            res(n,1) = res(ii,1);
            res(n,3) = adj_table(ii,1);
            res(n,4) = res(n-1,5);
            res(n,5) = res(n,4) + ( res(n,7) - res(n,6) );
            if ( res(n,4) > lim_range )
                res(n,4) = 0;
            end
            if ( res(n,5) > lim_range )
                res(n,5) = 0;
            end
            n = n + 1;
        end
    end
end

% --------------------------------------------------
% FUNCTIONS:
% --------------------------------------------------

    function d = fnd_clst( data, t_lim )

    %all_abs_mean = int64(mean(abs(data)));
    all_abs_std = int64(std(abs(data)));

    all_abs_median = int64(median(abs(data))); % need to check for mode!

    n_clst = all_abs_std / all_abs_median;

    if (n_clst == 0)
        d{1}(:,1) = data;
        d{1}(:,2) = [1:numel(data)];
        disp('The number of detected clusters is 0.');
        return;
    end

    T = clusterdata( data, n_clst*5); % Matlab: constructing agglomerative clusters from data using n_clst*5 clusters

    unq = unique(T);

    vars = zeros(1,2);

    l = 0;
    for i2 = 1:numel(unq) % for every detected cluster calculate skew std within it
        j2 = find( T == unq(i2) );
        i_abs_mean = int64( mean( abs(data(j2)) ) );
        %if ( i_abs_mean <= all_abs_median )
        if ( i_abs_mean <= t_lim )
            l = l + 1;
            vars(l,1) = std( data(j2) ); % std of skews within the cluster indexed as i2
            vars(l,2) = unq(i2); % cluster index
        end
    end
    
    % we expecting std = 0 if only 1 value is in the group;
    % should this be additionally addressed ?
    r_non0 = find( vars(:,1) ~= 0.0 );

    [ ~, r ] = min( vars(r_non0,1) ); % among non-zero skew std select the one with lowest value
    
    % determine the two cluster groups
    icluster0 = vars( r_non0(r),2 ); % the lowest-std cluster index
    icluster1 = vars( vars(:,1) == 0.0, 2 ); % select whose clusters which have std = 0 (single value in a cluster)

    if ( ~isempty(icluster0) ) % if there is no lowest-std cluster
        icluster = zeros( numel(icluster1)+1,1 );
        icluster(1,1) = icluster0;
        icluster(2:size(icluster),1) = icluster1(:,1);
    else
        icluster = zeros( numel(icluster1),1 );
        icluster(1:size(icluster),1) = icluster1(:,1);
    end


    d = {};
    for i2 = 1:size(icluster) % rearrange data according to clustring
        r = find( T == icluster(i2) );
        d{i2}(:,1) = data(r);
        d{i2}(:,2) = r;
    end

    end

end
