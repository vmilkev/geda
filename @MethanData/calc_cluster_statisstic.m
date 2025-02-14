function s = calc_cluster_statisstic( this, data, n_clst )
    % c = this.cluster_data( data, n_clst, 1 );
    % unq = unique(c); % number of unique clusters
    % if unq ~= n_clst
    %     warning("The number of unique clusters is not equal to the requested number of clusters: unq ~= n_clst.");
    % end
    % opt_vars = zeros(n_clst,2);
    % for i = 1:numel(unq)
    %     j = find( c == unq(i) );
    %     opt_vars(i,1) = max(data(j));
    %     opt_vars(i,2) = ( max(data(j)) - mean(data(j)) )/std(data(j));
    % end
    % [~,ind] = max( opt_vars(:,1) );
    % s = opt_vars(ind,2);

    B = sort(data);
    bin_size = floor(size(B,1)/n_clst);
    if bin_size == 0
        warning("The data size is too small: bin_size == 0!");
        s = 0;
        return;
    end
    d = B(end-bin_size:end,1);
    s = ( max(d) - mean(d) )/std(d);

    % up_limit = mean(data) + 3.0*std(data);
    % f = find(data(:,1) > up_limit);
    % d = data(f,1);
    % s = ( max(d) - mean(d) )/std(d);
end