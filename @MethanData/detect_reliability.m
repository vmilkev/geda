function [ reliability_and_skew, optimization_data, n_clusters ] = detect_reliability( this, time_skew, use_matl_opt )

n_clusters = 0; % optimal number of clusters
max_clusters = size(time_skew,1); % num clusters will be tested; this depends on data size

optimization_data = zeros(max_clusters,4); % storage for additional (usefull) info
clusters_matrix = zeros(size(time_skew,1), max_clusters); % storage for collecting clustering variants

for i = 1:max_clusters % prepare data for finding optimal num clusters
    clusters = cluster_data(time_skew, i); % cluster data using i num clusters    
    rlb_list = assign_reliable( time_skew, clusters ); % do binary classification using i num clusters
    clusters_matrix(:,i) = clusters; % collect clustering results
    optimization_data(i,1) = i; % num clusters
    optimization_data(i,3) = sum(rlb_list); % num reliable detections
end

optimization_data(:,2) = evaluate_clusters( time_skew, clusters_matrix, use_matl_opt ); % calculate clustering opt measure

switch use_matl_opt % apply specific filtering
    case 5
        optimization_data = optimization_data( optimization_data(:,2) > 1e-9, : );
    case 4
        optimization_data = filtering_data( optimization_data, 2 );
    case {1,3,2}
        optimization_data = filtering_data( optimization_data, 2 );
        optimization_data = optimization_data( optimization_data(:,2) < 1e15, : );
    otherwise
        warning('Cannot determine clustring criterion parameter!');
        return;
end

if isempty(optimization_data)
    warning('WARNING: Not enough data to detect reliability using clustering optimization!');
    reliability_and_skew = [];
    return;
end

optimization_data(1,4) = optimization_data(1,2)/optimization_data(1,2);
for i = 2:size(optimization_data,1) % calculate rel. change in clustering opt. measure
    optimization_data(i,4) = ( optimization_data(i,2) - optimization_data(i-1,2) )/optimization_data(i,2);
end

switch use_matl_opt % find optimal number of clusters
    case 5
        [~, ind] = min(optimization_data(:,4));
    case 4
        [~, ind] = min(optimization_data(:,4));
        [m_val, ind] = max(optimization_data(1:ind,2));
        f = find(optimization_data(1:ind,2) == m_val);
        ind = f(end);
    case {1,2,3}
        [m_val, ind] = max(optimization_data(2:end,4));
        if isempty(m_val)
            warning('WARNING: Not enough data to detect reliability using clustering optimization!');
            reliability_and_skew = [];
            return;
        end
        %ind = ind + 1; % because we apply max to the range 2:end
        f = find( abs( optimization_data(:,4) - m_val ) < 1e-2); % on some data gives better ressult
        ind = f(end);
    otherwise
        warning('Cannot determine clustring criterion parameter!');
        return;
end

n_clusters = optimization_data(ind,1); % optimal num clusters

clusters = cluster_data(time_skew, n_clusters); % cluster data using (optimal num) n_clusters
reliability_res = assign_reliable( time_skew, clusters ); % do reliability detection with optimized num. clusters

% Prepare results:
reliability_and_skew(:,1) = time_skew; % data
reliability_and_skew(:,2) = reliability_res; % reliability (binary classification) result

% NESTED FUNCTIONS:
    % ---------------------------------------------------------------------
    function res = evaluate_clusters( data, clst_matr, crit )
        rescaled_data = data;%rescale(data); % rescale data to avoid negative values
        switch crit
            case 1
                criter = 'CalinskiHarabasz'; % in average shows best accuracy (0.997)
            case 2
                criter = 'DaviesBouldin'; % still needs some work to locate optimum (at least two min often present)
            case 3
                criter = 'gap'; % (in log scale) is very similar to CalinskiHarabasz, but is terribly computationally demanding !!!
                eva_clst = evalclusters(rescaled_data,"kmeans","gap","KList",1:size(clst_matr,2));
                res = eva_clst.CriterionValues';
                return;
            case 4
                criter = 'silhouette'; % shows very good accuracy (0.941 aver.), but very inaccurate for some (one) data set!
            case 5
                criter = 'custom'; % best average accuracy 0.934
                [ res, ~ ] = evaluate_clusters_std( rescaled_data, clst_matr );
                return;
            otherwise
                warning("Cannot determine clustring criterion parameter!");
                return;
        end
        eva_clst = evalclusters(rescaled_data,clst_matr,criter);
        res = (eva_clst.CriterionValues');
    end
    % ---------------------------------------------------------------------
    function [ opt_measure, opt_vars ] = evaluate_clusters_std( data, clst_matr )        
        num_evaluat = size(clst_matr,2); % num clusters tested
        opt_measure = zeros( num_evaluat, 1 ); % cluster optimality measure, based on rel std        
        for j = 1:num_evaluat % loop over tested variants (diff num clusters tested)
            c = clst_matr(:,j);
            unq = unique(c); % number of unique clusters
            opt_vars = zeros(1,3);    
            l = 0;
            for i2 = 1:numel(unq)
                j2 = find( c == unq(i2) );
                if ( numel(j2) < 2 ) % initial filtering by excluding singe-valued clusters
                    continue;
                end
                i_rel_std = 0.0;
                if mean( data(j2) ) ~= 0.0
                    i_rel_std = std( data(j2) )/mean( data(j2) );
                end
                l = l + 1;
                opt_vars(l,1) = abs(i_rel_std); % relative std within the cluster
                opt_vars(l,2) = unq(i2); % cluster ID
                opt_vars(l,3) = numel(j2); % num records in cluster
            end    
            within_clst_dist = sum( opt_vars(:,1)./(opt_vars(:,3)) ); % weited sum of rel.std over all accepted clusters            
            opt_measure(j,1) = within_clst_dist;
        end
    end
    % ---------------------------------------------------------------------
    function [ classifcat_result ] = assign_reliable( data, clusters )
        [ ~, vars ] = evaluate_clusters_std( data, clusters ); % just getting vars matrix
        classifcat_result = zeros(size(data,1),1);
        for i2 = 1:size(vars,1) % do binary classification
            r = find( clusters == vars(i2,2) );
            classifcat_result(r,1) = 1; % TRUE only for those records which belong to clusters in vars
        end
    end
    % ---------------------------------------------------------------------
    function out_data = filtering_data( in_data, col )
        out_data = in_data;
        mask_val = isnan(out_data(:,col));
        out_data = out_data(~mask_val,:);
        mask_val = isinf(out_data(:,col));
        out_data = out_data(~mask_val,:);     
    end
    % ---------------------------------------------------------------------
    function c = cluster_data(data, n_clst)
        opts = statset('Display','off','MaxIter',500);
        c = kmeans(data,n_clst,'Options',opts,'emptyaction','singleton','replicate',5);
    end
    % ---------------------------------------------------------------------

end % end of the parent function
