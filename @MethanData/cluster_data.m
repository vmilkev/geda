function c = cluster_data( this, data, n_clst, rep )
opts = statset('Display','off','MaxIter',500);
c = kmeans(data,n_clst,'Options',opts,'emptyaction','singleton','replicate',rep);
end