function s = calc_cluster_statisstic( this, data, n_clst )
    B = sort(data);
    bin_size = floor(size(B,1)/n_clst);
    if bin_size == 0
        warning("The data size is too small: bin_size == 0!");
        s = 0;
        return;
    end
    d = B(end-bin_size:end,1);
    s = ( max(d) - mean(d) )/std(d);
end