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

% figure(1), plot(mesr_num{1, 1}(:,1),mesr_num{1, 1}(:,2),'o');
% figure(2), plot(mesr_num{1, 1}(:,1),mesr_num{1, 1}(:,3),'o');
% figure(3), plot(mesr_num{1, 1}(:,1),mesr_num{1, 1}(:,4),'o');
% figure(4), plot(mesr_num{1, 1}(:,1),mesr_num{1, 1}(:,5),'o');
% figure(5), plot(mesr_num{1, 1}(:,1),mesr_num{1, 1}(:,6),'o');
%% Reporting

% rp_file = "report_device_";
% rp_file1 = strcat(rp_file,num2str(m.device),".txt");
% rp_file2 = strcat(rp_file,num2str(m.device),".ps");
% 
% m.set_report_param( rp_file1, rp_file2 );

m.make_report("dat", "Time of reading all AMS data, seconds: ", floor(timeElapsed_r1) );
m.make_report("dat", "-------------------------------------------------------------------------------------------", []);
m.make_report("dat", " ", []);
m.make_report("txt", "AMS DATA", lely_num, lely_map);
m.make_report("dat", "Time of reading sniffer data, seconds: ", floor(timeElapsed_m1) );
m.make_report("dat", "Reading of robot: ", m.get_mesr_robot());
m.make_report("dat", "-------------------------------------------------------------------------------------------", []);
m.make_report("dat", " ", []);
m.make_report("txt", "SNIFFER DATA", mesr_num, mesr_map);
%m.make_report("fig", [], mesr_num, mesr_map, lely_num, lely_map);
%m.make_report("fig", [], mesr_num);

%% Clear some data
% we don't need these variables anymore
clear mesr_num mesr_map lely_num lely_map

%% Transform numeric robot record to Time Series data

% we've red and, hence, conside a specific robot:

if ( isempty(samplingT) || isnan(samplingT) ) % timing interval
    deltaT = 1.0;
else
    deltaT = samplingT; % timing interval
end

tic;
lely_ts = m.get_lely_ts( m.get_mesr_robot(), deltaT );
timeElapsed_r2 = toc;

%debug
% tdiff2 = lely_ts{1, 1}(2:end,1)-lely_ts{1, 1}(1:end-1,1);
% figure, plot(tdiff2*(24 * 60 * 60),'o' );

%% Adjust Merurements TS to deltaT value used in Robot TS

% since we have red only a specific robot mesurements:
tic;
[msr_ts,msr_ts_noscale] = m.get_mesr_ts(deltaT);
timeElapsed_m2 = toc;

% figure(6), plot(msr_ts_noscale{1, 1}(:,1),msr_ts_noscale{1, 1}(:,2), 'o');
% figure(7), plot(msr_ts_noscale{1, 1}(:,1),msr_ts_noscale{1, 1}(:,3), 'o');
% figure(8), plot(msr_ts_noscale{1, 1}(:,1),msr_ts_noscale{1, 1}(:,4), 'o');
% figure(9), plot(msr_ts_noscale{1, 1}(:,1),msr_ts_noscale{1, 1}(:,5), 'o');
% figure(10), plot(msr_ts_noscale{1, 1}(:,1),msr_ts_noscale{1, 1}(:,6), 'o');
%% Reporting

m.make_report( "dat", "Sampling interval used to transform measured records to time series, seconds: ", deltaT );
m.make_report( "dat", "Elapsed time of transforming AMS records to time series, seconds:             ", floor(timeElapsed_r2) );
m.make_report( "dat", "Elapsed time of transforming sniffer records to time series, seconds:         ", floor(timeElapsed_m2) );
m.make_report("dat", "-------------------------------------------------------------------------------------------", []);
%m.make_report("fig", [], msr_ts, lely_ts);

%% Filtering step

rfirst = 1;

A = lely_ts{1,1}(1:end,2);

% A2 = zeros(size(A));
% Acell = {};
% 
% if is_ssa % do SSA on A
% 
%     w_length = 4000;
%     n_waves = floor( numel(A)/w_length ); % number of waves (sub-signals)
% 
%     r1 = zeros(n_waves);
%     r2 = zeros(n_waves);
%     r1(1) = 1;
%     r2(1) = w_length;
%     for ii = 2:n_waves
%         r1(ii) = r1(ii-1) + w_length;
%         r2(ii) = r2(ii-1) + w_length;
%     end
%     r2(n_waves) = numel(A);
% 
%     parfor ii = 1:n_waves
%         [A0,~] = m.ssa2(A(r1(ii):r2(ii)),ssa_wnd,ssa_eigval);
%         Acell{ii} = A0;
%     end
% 
%     for ii = 1:n_waves
%         A2( r1(ii):r2(ii) ) = Acell{ii};
%     end
% end

if ~is_ssa
    A = (A - 0.5)*2.0; % scalling to [-1,1]
end

test_statistic = zeros(numel(msr_ts),1);

res_elapsed = {};

skew = zeros(numel(msr_ts),1);

adj_table = zeros(numel(msr_ts),5);

if ( isempty(signallength) || isnan(signallength) ) % number of signals to detect
    sig_length = 5.0; % is an assumed minimal duration of a signal, [hours]
else
    sig_length = signallength;
end

aligned_data_file = "aligned_init_sniffer_";
aligned_data_file = strcat(aligned_data_file,num2str(m.device),".txt");

aligned_data_file2 = "aligned_adj_sniffer_";
aligned_data_file2 = strcat(aligned_data_file2,num2str(m.device),".txt");

aligned_data_file3 = "aligned_rlb_sniffer_";
aligned_data_file3 = strcat(aligned_data_file3,num2str(m.device),".txt");

if exist(aligned_data_file, 'file')
    delete(aligned_data_file);
end

if exist(aligned_data_file2, 'file')
    delete(aligned_data_file2);
end

if exist(aligned_data_file3, 'file')
    delete(aligned_data_file3);
end

varNames = {'time','id','gas','signal'};

all_signals = 1;

for j = 1:numel(msr_ts) % for each Sniffer data file

    if ( size(msr_ts{j,1},1) <= 1  )
        m.make_report("dat", "WARNING! There is only 1 or no records present in the sub-dataset no.:", j);
        m.make_report("dat", "-------------------------------------------------------------------------------------------", []);
        continue;
    end

    snif_length = etime( datevec(msr_ts{j,1}(end,1)),datevec(msr_ts{j,1}(1,1)) )/( 60*60 ); % duration of sniffer records, [hours]

    if ( snif_length < 1 )
        m.make_report("dat", "WARNING! The signal length (duration) is less than 1 hour in the sub-dataset no.:", j);
        m.make_report("dat", "         This sub-dataset will not be processed.", []);
        m.make_report("dat", "-------------------------------------------------------------------------------------------", []);
        continue;
    end
    
    waves = floor( snif_length/sig_length ); % number of waves (signals)

    if ( waves == 0 )
        waves = 1;
        m.make_report("dat", "WARNING! The signal length (duration) provided with $SIGLENGTH parameter is too high for the processing data!", []);
        m.make_report("dat", "         Entire sniffer records will be processed without splitting!", []);
        m.make_report("dat", "-------------------------------------------------------------------------------------------", []);
    end

%     vars = string([1:size(msr_ts_noscale{j,1},2)]');
%     for i_vars =1:size(msr_ts_noscale{j,1},2)
%         vars(i_vars) = strcat('var_',vars(i_vars));
%     end
        
    wlength = floor(numel(msr_ts{j,1}(:,2))/waves);

    tic;
    for k = 1:waves
        
        first = 1+(k-1)*wlength;
        last = k*wlength;

        if ~is_ssa
            B = msr_ts{j,1}(first:last,2);
        else
            B = msr_ts_noscale{j,1}(first:last,3); % No scalling
            [B,~] = m.ssa2(B,ssa_wnd,ssa_eigval);
        end
            
        range = numel(A(:,1))-numel(B(:,1))+1;
        
        res = zeros( range,1 );
        res2 = zeros( range,1 );

        el_B = numel(B(:,1));
    
        parfor i = 1:range
            % not really convolution, just one step of it, at the expected max overlap. Calculating only maximums of convolution output            
            res(i,1) = dot( A(i:i+el_B-1),B )/(norm(A(i:i+el_B-1))*norm(B)); % Test statistic
        end

        [~,r] = max( res ); % version 2

        test_statistic(j,k) = res(r,1); % test statistic
            
        lely_time = datestr( lely_ts{1,1}( r,1 ) );
        
        snifer_time = datestr( msr_ts{j,1}(first,1) );

        if ( size(snifer_time,2) < 15 )
            snifer_time = strcat(snifer_time, ' 00:00:00');
            disp('alignTS(): found short sniffer time.');
        end

        ams_dates{j,k} = lely_time;
        snifer_dates{j,k} = snifer_time;
        
        t = string( datestr( lely_ts{1,1}( r:r+el_B-1,1 ) ) );
        id = lely_ts{1,1}( r:r+el_B-1,3 );
        gas = msr_ts_noscale{j,1}(first:last,2:end);
        
        % signal number
        signal_no = zeros( numel(B(:,1)),1 );
        signal_no(:,1) = all_signals;

        write_header = true;
        if ( all_signals > 1 )
            write_header = false;
        end

        if (deltaT >= m.deltaT_original)
            refine = ~mod( uint64( (1:(last-first)+1) )', 1 );
        else
            refine = ~mod( uint64( (1:(last-first)+1) )', floor(m.deltaT_original/deltaT) );
        end

        aligned = table( t(refine,:),id(refine,:),gas(refine,:),signal_no(refine,:),'VariableNames',varNames );

        writetable(aligned,aligned_data_file,'WriteMode','Append',...
                    'WriteVariableNames',write_header,'Delimiter',';','QuoteStrings',true);
                
        skew(j,k) = seconds(diff(datetime([lely_time;snifer_time])));

        % save information of each alignment,
        % will be used for correction/adjustment
        % for an incorrectly alighned data
        adj_table(all_signals,1) = skew(j,k);
        adj_table(all_signals,2) = r;
        adj_table(all_signals,3) = r+el_B-1;
        adj_table(all_signals,4) = first;
        adj_table(all_signals,5) = last;
        adj_table(all_signals,6) = j;

        all_signals = all_signals + 1;
    
    end
    
    res_elapsed{j,1} = toc; %this variant for the paper: toc/waves;
    
    m.make_report("dat", "Completed. Processing of the sniffer data set no.:                                ", j);
    m.make_report("dat", "           Elapsed time of filtering all signals sourced from this file, seconds: ", floor(res_elapsed{j,1}) );
    m.make_report("dat", "           Number of skew estimates (and, hence, filtered signals):               ", waves);
    m.make_report("dat", "           Duration (length) of entire sniffer records, hours:                    ", floor(snif_length) );
    m.make_report("dat", "           Duration (length) of each signal, hours:                               ", sig_length);
    m.make_report("dat", "-------------------------------------------------------------------------------------------", []);

end

stats{1} = skew;
stats{2} = test_statistic;

% stats{2} = signal;
% stats{3} = signalstd;
% stats{4} = test_statistic;
% stats{5} = signalqual;
% stats{6} = signalqual2;

%% Reporting
m.make_report("stats", "Skews and MF signals:", stats);
m.make_report("skew", "Sniffer starting dates (for each signal):", snifer_dates, ams_dates);

%% CORRECTING (ADJUSTING) ALIGNMENT OF UNRELIABLE DATA

m.make_report("dat", "Correcting (re-adjusting) the unreliable data.", []);
m.make_report("dat", "-------------------------------------------------------------------------------------------", []);

res = m.adjust_skews( adj_table, numel(A(:,1)), sig_length*60*60 );

if ( isempty(res) )
    m.make_report("dat", " ", []);
    m.make_report("dat", "WARNING! All processed data found to be UNRELIABLE. There is no successful alignement.", []);
    m.make_report("dat", "         Correction (re-adlustment) is not possible. The data should be revised manually.", []);
    m.make_report("dat", " ", []);
    return;
end

if ( all(res(:,2)) )
    m.make_report("dat", " ", []);
    m.make_report("dat", "Completed. All processed data found to be RELIABLE.", []);
    m.make_report("dat", "           No further correction (re-adlustment) is needed.", []);
    m.make_report("dat", " ", []);
    return;
end

for i = 1:size(res,1)
    
    a1 = res(i,4);
    a2 = res(i,5);
    b1 = res(i,6);
    b2 = res(i,7);
    j = res(i,8);

    if ( a1 == 0 || a2 == 0 || b1 == 0 || b2 == 0 )
        continue;
    end

    t = string( datestr( lely_ts{1,1}( a1:a2,1 ) ) );
    id = lely_ts{1,1}( a1:a2,3 );
%     co2 = msr_ts_noscale{j,1}(b1:b2,3);
%     ch4 = msr_ts_noscale{j,1}(b1:b2,2);
    gas = msr_ts_noscale{j,1}(b1:b2,2:end);

    % signal number
    signal_no = zeros( a2-a1+1,1 );
    signal_no(:,1) = i;

    %aligned = table(t,id,co2,ch4,signal_no);
    %aligned = table(t,id,gas,signal_no);

    if (deltaT >= m.deltaT_original)
        refine = ~mod( uint64( (1:(b2-b1)+1) )', 1 );
    else
        refine = ~mod( uint64( (1:(b2-b1)+1) )', floor(m.deltaT_original/deltaT) );
    end

    aligned = table( t(refine,:),id(refine,:),gas(refine,:),signal_no(refine,:) );

    write_header = true;
    if ( i > 1 )
        write_header = false;
    end

    writetable(aligned,aligned_data_file2,'WriteMode','Append',...
        'WriteVariableNames',write_header,'Delimiter',';','QuoteStrings',true);

    if (res(i,2))
        writetable(aligned,aligned_data_file3,'WriteMode','Append',...
            'WriteVariableNames',write_header,'Delimiter',';','QuoteStrings',true);
    end

    lely_time = datestr( lely_ts{1,1}( a1,1 ) );
    snifer_time = datestr( msr_ts{j,1}(b1,1) );

    if ( size(snifer_time,2) < 15 )
        snifer_time = strcat(snifer_time, ' 00:00:00');
    end

    ams_dates2{1,i} = lely_time;
    snifer_dates2{1,i} = snifer_time;

end

m.make_report("dat", "Completed. All processed data successfully re-adjusted.", []);
m.make_report("dat", "-------------------------------------------------------------------------------------------", []);
m.make_report("dat", " ", []);
m.make_report("stats", "Skews and MF signals:", res, adj_table);
m.make_report("skew", "Re-adjusted sniffer starting dates (for each signal):", snifer_dates2, ams_dates2);

delete(m);

end