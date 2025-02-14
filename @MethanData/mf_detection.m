function [res_stats, adj_tab, skew_res, test_stat1, test_stat2, ams_dt, snf_dt] = mf_detection(this, snf_data, ams_data, ams_arr, sig_len, use_ssa, report_msg)

test_stat1 = cell(numel(snf_data),1);
test_stat2 = cell(numel(snf_data),1);

skew_res = cell(numel(snf_data),1);

adj_tab = zeros(numel(snf_data),7);

all_signals = 1;

res_stats = zeros(numel(snf_data), 5);

for j = 1:numel(snf_data) % for each Sniffer data file or data set

    if ( size(snf_data{j,1},1) <= 1  )
        if report_msg
            this.make_report("dat", "----------------------------------------------------------------------", []);
            this.make_report("dat", "WARNING! There is only 1 or no records present in the signal no.:", j);
            this.make_report("dat", "----------------------------------------------------------------------", []);
        end
        continue;
    end

    snif_length = etime( datevec(snf_data{j,1}(end,1)),datevec(snf_data{j,1}(1,1)) )/( 60*60 ); % duration of sniffer records, [hours]

    if ( snif_length < 1 )
        if report_msg
            this.make_report("dat", "--------------------------------------------------------------------------------------------", []);
            this.make_report("dat", "WARNING! The signal length is less than 1 hour, hence will not be processed; the singal no.:", j);
            this.make_report("dat", "--------------------------------------------------------------------------------------------", []);
        end
        continue;
    end

    waves = floor( snif_length/sig_len ); % number of waves (signals)

    if ( waves == 0 )
        waves = 1;
        % this.make_report("dat", "WARNING! The signal length (duration) provided with $SIGLENGTH parameter is too high for the processing data!", []);
        % this.make_report("dat", "         Entire sniffer records will be processed without splitting!", []);
        % this.make_report("dat", "-------------------------------------------------------------------------------------------", []);
    end

    wlength = floor(numel(snf_data{j,1}(:,2))/waves);

    tic;
    for k = 1:waves % for each signal of length (duration) sig_len

        first = 1+(k-1)*wlength;
        last = k*wlength;

        if ~use_ssa
            B = snf_data{j,1}(first:last,2);
        else
            B = snf_data{j,1}(first:last,2);
            [B,~] = this.ssa2(B,ssa_wnd,ssa_eigval);
        end

        range = numel(ams_arr(:,1))-numel(B(:,1))+1;

        res = zeros( range,1 );

        el_B = numel(B(:,1));

        nrm_B = norm(B);

        parfor i = 1:range
            res(i,1) = dot( ams_arr(i:i+el_B-1),B )/( norm(ams_arr(i:i+el_B-1))*nrm_B ); % collecting test statistic
        end

        [~,r] = max( res );

        test_stat1{j}(k) = res(r,1); % test statistic
        test_stat2{j}(k) = this.calc_cluster_statisstic( res, 100 );

        lely_time = datetime( ams_data{1,1}( r,1 ),'ConvertFrom','datenum','Format','dd-MM-yyy HH:mm:ss' );
        snifer_time = datetime( snf_data{j,1}(first,1),'ConvertFrom','datenum','Format','dd-MM-yyy HH:mm:ss' );

        if ( size(char(snifer_time),2) < 15 ) % the case when there is no HH:mm:ss in snifer_time
            snifer_time = strcat( snifer_time, ' 00:00:00');
            if report_msg
                this.make_report("dat", "-------------------------------------------------------------------------------------", []);
                this.make_report("dat", "WARNING: found short sniffer time; the case when there is no HH:mm:ss in snifer_time.", []);
                this.make_report("dat", "-------------------------------------------------------------------------------------", []);
            end
        end

        ams_dt{j,k} = lely_time;
        snf_dt{j,k} = snifer_time;

        skew_val = round( seconds(diff(datetime([lely_time;snifer_time]))) );
        skew_res{j}(k) = skew_val;

        % save information of each alignment, will be used for reliability detection
        adj_tab(all_signals,1) = skew_val; % calculated skew_res
        adj_tab(all_signals,2) = r; % starting index of the signal k in the AMS TS
        adj_tab(all_signals,3) = r+el_B-1; % ending index of the signal k in the AMS TS
        adj_tab(all_signals,4) = first; % starting index of the signal k in the sniffer TS
        adj_tab(all_signals,5) = last; % ending index of the signal k
        adj_tab(all_signals,6) = j; % index of the processued dataset (often a file)

        all_signals = all_signals + 1;

    end

    elps_time = toc;

    res_stats(j,1) = j; % data set number
    res_stats(j,2) = floor(snif_length); % duration of the dataset, hours
    res_stats(j,3) = waves; % num of skew_res estimates
    res_stats(j,4) = floor( snif_length/waves ); % duration of each estimate
    res_stats(j,5) = elps_time; % elapsed time of each estimate

end

end