function alignTS( fparam )

% fparam - parameter file;

%% Create MethanData class object

m = MethanData();

%% Open report (log) file

m.set_report_param( strcat("log_geda_",m.date_now,".txt") );

%% Set GEDA parameters

[ samplingT, signallength ] = m.set_geda_param(fparam);

%% extract data in the numeric/table format

tic;
[lely_num, ~] = m.get_lely_data( "numeric" );
[lely_map, ~] = m.get_lely_data( "table" );
timeElapsed_r1 = toc;

tic;
mesr_num = m.get_mesr_data( "numeric" );
mesr_map = m.get_mesr_data( "table" );
timeElapsed_m1 = toc;

%% Checking if the Sniff data is within the range of AMS data

for i = 1:size(mesr_num,1)
    if ( lely_num{1, 1}(end,2)-lely_num{1, 1}(1,2) ) < ( mesr_num{i, 1}(end,1)-mesr_num{i, 1}(1,1) )
       m.make_report("datint", "WARNING: AMS data range is smaller than SNIFF data range. The data set no.:", i);
    end
    if lely_num{1, 1}(1,2) >= mesr_num{i, 1}(end,1)
        m.make_report("datint", "WARNING: AMS data range starts later than the SNIFF data range ends. The data set no.:", i);
    end
    if lely_num{1, 1}(end,2) <= mesr_num{i, 1}(1,1)
        m.make_report("datint", "WARNING: AMS data range ends earlier than the SNIFF data range starts. The data set no.:", i);
    end
end

%% Reporting

m.make_report("dat", "Time of reading all AMS data, seconds: ", timeElapsed_r1 );
m.make_report("dat", " ", []);
m.make_report("txt", "AMS DATA", lely_num, lely_map);
m.make_report("dat", "Time of reading sniffer data, seconds: ", timeElapsed_m1 );
m.make_report("datint", "Reading of robot: ", m.get_mesr_robot());
m.make_report("dat", " ", []);
m.make_report("txt", "SNIFFER DATA", mesr_num, mesr_map);

%% Clear some data

clear mesr_num mesr_map lely_num lely_map % we don't need these variables anymore

%% Set sampling frequency

if ( ~isempty(samplingT) && ~isnan(samplingT) ) 
    m.set_sampling_freq(samplingT); % set sampling frequency if provided by the user
end

%% Transform numeric robot record to Time Series data

% msr_ts - rescaled and resampled TS with two cols: 1 - serial time; 2 - CO2
% snf_ts_fnames - file names of not sacled but resampled TS data with all
% cols as in the original sniffers data, though, the col 1 - is serial time

[ msr_ts, snf_ts_fnames ] = m.get_mesr_ts();

% lely_ts - not scaled but resampled AMS data with three cols: 1 - serial
% time; 2 - value {0,1}; 3 - IDs.
% ams_fnames - binary file names for resampled at 1 sec frequency ams data

[ lely_ts, ams_fnames ] = m.get_lely_ts( m.get_mesr_robot() );


%% Prepare some data

A = lely_ts{1,1}(1:end,2);
A = (A - 0.5)*2.0; % scalling to [-1,1]

if m.is_ssa % do SSA on A
    tic;
    m.make_report("dat", "Starting SSA on AMS data...", []);
    max_length = 4000;
    A = m.apply_ssa_ams(A, min(numel(A),max_length), m.ssa_wnd, m.ssa_eigval);
    elps_t = toc;
    m.make_report("dat", "SSA on AMS data is completed. Elapsed time: ", elps_t);
end

%% Set the sig_length parameter or calculate an optimal value of sig_length

if ( signallength == 0 ) % if no parameter is provided, use the default one
    sig_length = 8;
elseif signallength < 0 % in this case do optimization for the sig_length parameter
    max_sig_length = 8.0; % is an assumed minimal duration of a signal, [hours]
    snif_len = zeros(numel(msr_ts),1);
    for i = 1:numel(msr_ts) % find lengths of sniffer records
        snif_len(i) = etime( datevec(msr_ts{i,1}(end,1)),datevec(msr_ts{i,1}(1,1)) )/( 60*60 ); % duration of sniffer records, [hours]
    end
    min_snif_len = mode(snif_len); % find length of smallest data set
    opt_res = zeros(min( floor(min_snif_len), max_sig_length),7); % reliability detection results for each tested sig. length
    l = 1;
    for i_len = 2:min(min_snif_len, max_sig_length) % suppose we use as the upper limit of signal duration is either 5 hours or length of smallest data set
        tic;
        [~,adj_table,~,~,~,~,~] = m.mf_detection(msr_ts, lely_ts, A, i_len, m.is_ssa, false);
        [ rlb_res, ~, ~ ] = m.detect_reliability(adj_table(:,1), 1);
        if isempty(rlb_res)
            rlb_res = m.detect_reliability2( adj_table, i_len*60*60 );
        end
        [rlb_confidence, corected_rlb] = m.get_rlb_accuracy( rlb_res );
        elp_time = toc;
        f = find(corected_rlb(:,1) == 1);
        opt_res(i_len-l,1) = i_len; % sig. length tested
        opt_res(i_len-l,2) = 100*sum(rlb_res(:,2))/numel(rlb_res(:,2)); % prpportion of reliabile data
        opt_res(i_len-l,3) = 100*rlb_confidence; % confidence on reliabile data
        opt_res(i_len-l,4) = mean(adj_table(f,1)); % mean skew
        opt_res(i_len-l,5) = std(adj_table(f,1)); % std skew
        opt_res(i_len-l,6) = elp_time; % elapsed time
        % opt_res(i_len-l,7) = 0; - this is dummy cols needed for reporting

        if opt_res(i_len-l,2) == 100
            break;
        end
    end
    [~,ind] = max(opt_res(:,2));
    sig_length = opt_res(ind,1);

    m.make_report("stats", "Summary of signal length optimization:", opt_res);
else
    sig_length = signallength;
end

m.make_report("datint", "Maximal signal length used for MF detection, hours:", sig_length);
m.make_report("dat", " ", []);

%% Clear unnecessary data

clear opt_res adj_table rlb_res snif_len

%% MF detection

[processing_stats,...
 adj_table,...
 skew,...
 test_statistic,...
 test_statistic_clst,...
 ams_dates,...
 snifer_dates] = m.mf_detection(msr_ts, lely_ts, A, sig_length, m.is_ssa, true);

%% Reliability estimation of performed detection

[ rlb_res, ~, ~ ] = m.detect_reliability(adj_table(:,1), 1); % use CalinskiHarabasz criterion for clustering optimization

if isempty(rlb_res)
    rlb_res = m.detect_reliability2( adj_table, sig_length*60*60 );
end

[accuracy, corr_rlb] = m.get_rlb_accuracy( rlb_res );

%% Collect summary

if accuracy < 1.0
    adj_table(:,7) = corr_rlb(:,1);
else
    adj_table(:,7) = rlb_res(:,2);
end

rlb_evaluation = {size(skew,1)};
l = 0;
for i = 1:size(skew,1)
    rlb_evaluation{i,1} = adj_table( l+1:l+size(skew{i, 1},2 ), 7)';
    l = l + size(skew{i, 1},2);
end

stats{1,1} = skew;
stats{2,1} = test_statistic;
stats{3,1} = test_statistic_clst;
stats{4,1} = rlb_evaluation;

%% Reporting

m.make_report("stats", "Summary of MF detection:", processing_stats);
m.make_report("dat", "Data quality: proportion of reliable data,      %:", 100*sum(adj_table(:,7))/numel(adj_table(:,7)) );
m.make_report("dat", "              confidence of quality estimation, %:", accuracy*100 );
m.make_report("dat", " ", []);
m.make_report("stats", "Skews and MF signals:", stats);
m.make_report("skew", "Sniffer starting dates (for each signal):", snifer_dates, ams_dates);

%% Clear unnecessary data

clear A msr_ts lely_ts res stats snifer_dates ams_dates

%% Writing alignment results

tic;
if (sum(adj_table(:,7))/numel(adj_table(:,7))) > 0.0 % if the reliable data found
    m.write_alignment_results( adj_table, ams_fnames, snf_ts_fnames );
else
    m.make_report("dat", "WARNING: No reliable data found. There is nothing to write into files.", []);
end
elps_time = toc;

m.make_report("dat", "Elapsed time of writing the aligned data into files, sec:", elps_time);

%% Calculating and Writing emission traits

tic;
if (sum(adj_table(:,7))/numel(adj_table(:,7))) > 0.0 % if the reliable data found
    m.write_traits( adj_table, ams_fnames, snf_ts_fnames, sig_length );
else
    m.make_report("dat", "WARNING: No reliable data found. Emission traits will not be provided.", []);
end
elps_time = toc;

m.make_report("dat", "Elapsed time of calculating traits and writing into files, sec:", elps_time);
m.make_report("dat", " ", []);
print_type = "The produced trait is the gas emission intensity (based on every AMS visit calculations), [ppm/min].";
if m.trait_type == 2
    print_type = "The produced trait is the mean gas emission per visit (averaged over number of observations in an AMS visit), [ppm].";
end
m.make_report("dat", print_type, []);

%% Finalize: delete .bin files and close the log file

m.fclean(ams_fnames);
m.fclean(snf_ts_fnames);

m.make_report( "dat", " ", []);
m.make_report( "dat", strcat("GEDA LOG, completed: ",string(datetime()) ), []);

m.on_exit(); % close log file

return;

end