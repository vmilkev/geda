function d = detect_reliable( filename )

skew = read_skew(filename);

d = fnd_clst( skew, 1.3 );

% NESTED FUNCTIONS:

    function d = fnd_clst( data, clst_thr )

        all_abs_mean = int64(mean(sqrt(data.*data)));
        all_abs_median = int64(median(sqrt(data.*data))); % need to check for mode!
        n_clst = ceil(all_abs_mean / all_abs_median);

        if (n_clst == 0)
            d = zeros(size(data,1),2);
            d(:,1) = data;
            return;
        end
        
        % Pairwise distance between pairs of data

        %Y = pdist(data,'euclidean');
        %Y = pdist(data,'seuclidean');
        Y = pdist(data,'mahalanobis');
        % Y4 = pdist(data,'cityblock');
        % Y5 = pdist(data,'minkowski');
        % Y6 = pdist(data,'chebychev');
        
        % Z2 = linkage(Y,'centroid'); % Centroid distance (UPGMC), appropriate for Euclidean distances only
        % Z4 = linkage(Y,'median'); % Weighted center of mass distance (WPGMC), appropriate for Euclidean distances only
        % Z6 = linkage(Y,'ward'); % Inner squared distance (minimum variance algorithm), appropriate for Euclidean distances only
        
        Z = linkage(Y,'average'); % Unweighted average distance (UPGMA)
        Z3 = linkage(Y,'complete'); % Farthest distance
        Z5 = linkage(Y,'single'); % Shortest distance      
        Z7 = linkage(Y,'weighted'); % Weighted average distance (WPGMA)

        %dendrogram(Z);
        
        % T2 = cluster(Z2,'MaxClust',n_clst);
        % T4 = cluster(Z4,'MaxClust',n_clst);
        % T6 = cluster(Z6,'MaxClust',n_clst);
        
        T3 = cluster(Z3,'MaxClust',n_clst);
        T1 = cluster(Z,'MaxClust',n_clst);
        T5 = cluster(Z5,'MaxClust',n_clst);
        T7 = cluster(Z7,'MaxClust',n_clst);

        T = T1;
        
        unq = unique(T);

        vars = zeros(1,2);

        l = 0;
        for i2 = 1:numel(unq)
            j2 = find( T == unq(i2) );
            i_rel_std = std( rescale(data(j2)) + 1 )/mean( rescale(data(j2)) + 1 ); % rescaling because of possible negative numbers; +1 to avoid NAN if single0values clusters exist
            if ( i_rel_std <= 0.3 ) % initial filtering allowing no clusters found in case of unreliable data
                l = l + 1;
                vars(l,1) = i_rel_std; % relative std within the cluster
                vars(l,2) = unq(i2); % cluster ID
            end
        end

        r_non0 = find( vars(:,1) ~= 0.0 ); % exclude singe-value clusters
        [ ~, r ] = min( vars(r_non0,1) ); % cluster with min std within it
        std_threshold = vars(r,1) * clst_thr; % establish threshold to select clusters beow it
        r_thr = find( vars(:,1) <= std_threshold ); % find clusters below threshold

        d = zeros(size(data,1),2);
        d(:,1) = data;

        for i2 = 1:size(r_thr,1)
            r = find( T == vars(r_thr(i2),2) );
            d(r,2) = 1;
        end

    end

    function s = read_skew(filename, startRow, endRow)

        % Initialize variables.
        if nargin<=2
            startRow = 1;
            endRow = inf;
        end

        formatSpec = '%*10s%12f%[^\n\r]';

        fileID = fopen(filename,'r');

        if fileID == -1
            disp( strcat("Cannot open the file: ", filename) );
            return;
        end

        % Read columns of data according to the format.
        dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for block=2:length(startRow)
            frewind(fileID);
            dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            dataArray{1} = [dataArray{1};dataArrayBlock{1}];
        end

        fclose(fileID);

        s = [dataArray{1:end-1}];

    end

end % end of the parent function
