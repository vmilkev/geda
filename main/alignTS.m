function alignTS( fparam )

% fparam - parameter file;

%% Create MethanData class object

m = MethanData();

%% Set general GEDA parameters

[ samplingT, signallength, is_ssa, ssa_wnd, ssa_eigval ] = m.set_geda_param(fparam);

if isempty(is_ssa)
    is_ssa = false;
else
    if is_ssa == 0
       is_ssa = false;
    else
        is_ssa = true;
        if isempty(ssa_wnd)
            ssa_wnd = 200;
        end
        if isempty(ssa_eigval)
            ssa_eigval = 10;
        end
        if ssa_wnd < ssa_eigval
            disp('Parameter SSA1 cannot be lower than SSA2!');
            return;
        end
    end
end

%% Set parameters for reading data

m.set_lely_param(fparam);
m.set_mesr_param(fparam);

%% Set report

rp_file = "log_device_";
rp_file1 = strcat(rp_file,num2str(m.device),".txt");
rp_file2 = strcat(rp_file,num2str(m.device),".ps");

m.set_report_param( rp_file1, rp_file2 );

%% extract data in the numeric/table format

tic;
[lely_num, ~] = m.get_lely_data( "numeric" );
[lely_map, ~] = m.get_lely_data( "table" );
timeElapsed_r1 = toc;

tic;
mesr_num = m.get_mesr_data( "numeric" );
mesr_map = m.get_mesr_data( "table" );
timeElapsed_m1 = toc;

%% Reporting

m.make_report("dat", "Time of reading all AMS data, seconds: ", timeElapsed_r1 );
m.make_report("dat", " ", []);
m.make_report("txt", "AMS DATA", lely_num, lely_map);
m.make_report("dat", "Time of reading sniffer data, seconds: ", timeElapsed_m1 );
m.make_report("dat", "Reading of robot: ", m.get_mesr_robot());
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

if is_ssa % do SSA on A
    A2 = zeros(size(A));
    Acell = {};
    m.make_report("dat", "Starting SSA on AMS data...", []);
    A = (A - 0.5)*2.0; % scalling to [-1,1]
    w_length = 4000;
    n_waves = floor( numel(A)/w_length ); % number of waves (sub-signals)
    r1 = zeros(n_waves);
    r2 = zeros(n_waves);
    r1(1) = 1;
    r2(1) = w_length;
    for ii = 2:n_waves
        r1(ii) = r1(ii-1) + w_length;
        r2(ii) = r2(ii-1) + w_length;
    end
    r2(n_waves) = numel(A);

    parfor ii = 1:n_waves
        [A0,~] = m.ssa2(A(r1(ii):r2(ii)),ssa_wnd,ssa_eigval);
        Acell{ii} = A0;
    end
    for ii = 1:n_waves
        A2( r1(ii):r2(ii) ) = Acell{ii};
    end
    A = A2;
    m.make_report("dat", "SSA on AMS data is completed.", []);
else
    A = (A - 0.5)*2.0; % scalling to [-1,1]
end

%% Calculate optimal sig_length

if ( isempty(signallength) || isnan(signallength) ) % number of signals to detect
    max_sig_length = 6.0; % is an assumed minimal duration of a signal, [hours]
    snif_len = zeros(numel(msr_ts),1);
    for i = 1:numel(msr_ts) % find lengths of sniffer records
        snif_len(i) = etime( datevec(msr_ts{i,1}(end,1)),datevec(msr_ts{i,1}(1,1)) )/( 60*60 ); % duration of sniffer records, [hours]
    end
    min_snif_len = mode(snif_len); % find length of smallest data set
    opt_res = zeros(min( floor(min_snif_len), max_sig_length),7); % reliability detection results for each tested sig. length
    l = 1;
    for i_len = 2:min(min_snif_len, max_sig_length) % suppose we use as the upper limit of signal duration is either 5 hours or length of smallest data set
        tic;
        [~,adj_table,~,~,~,~,~] = m.mf_detection(msr_ts, lely_ts, A, i_len, is_ssa, false);
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

m.make_report("dat", "Maximal signal length used for MF detection, hours:", sig_length);
m.make_report("dat", " ", []);

clear opt_res adj_table rlb_res snif_len

%% MF detection

[processing_stats,...
 adj_table,...
 skew,...
 test_statistic,...
 test_statistic_clst,...
 ams_dates,...
 snifer_dates] = m.mf_detection(msr_ts, lely_ts, A, sig_length, is_ssa, true);

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
%m.make_report("dat", " ", []);

%% Clear unnecessary data

clear A msr_ts lely_ts res stats snifer_dates ams_dates

%% Writing alignment results

m.make_report("dat", "Writing the aligned data into files ...", []);

tic;
if (sum(adj_table(:,7))/numel(adj_table(:,7))) > 0.0 % if the reliable data found
    m.write_alignment_results( adj_table, ams_fnames, snf_ts_fnames );
end
elps_time = toc;

m.make_report("dat", "Completed. Elapsed time of writing the aligned data into files, sec:", elps_time);

%% delete .bin files

m.fclean(ams_fnames);
m.fclean(snf_ts_fnames);

delete(m); % explicitelly calling the class destructor; though, not really necessary

return;

end