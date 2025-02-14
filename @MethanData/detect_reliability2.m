function res = detect_reliability2( this, adj_table, lim_time )

% adj_table - alignment results matrix
% lim_time - emission signal length, sec 

s = adj_table(:,1); % estimated time skews

d = fnd_clst( s, lim_time );

res = zeros( size(s,1), 2 );

if ( ~any(d{1}) )
    disp('There are no clusters found in the data to be adjusted.');
    return;
end

for i = 1:size(d,2)
    j_correct = d{i}(:,2);
    res( j_correct, 1 ) = d{i}(:,1); % correct timing skews
    res( j_correct, 2 ) = 1; % flags for correct timing skews
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
        warning('The number of detected clusters is 0.');
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
