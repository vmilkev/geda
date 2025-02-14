function write_alignment_results( this, summary_table, ams_fnames, snf_fnames )

    date_now = regexprep(regexprep(string(datetime()), ' +', '_'), ':+', '-'); % day-time without sopaces and colomns in between
    aligned_data_file_rlb = strcat("rlb_",...
                                   date_now,...
                                   "_dev_",...
                                   num2str(this.device),... % robot (device) id
                                   ".geda");

    if exist(aligned_data_file_rlb, 'file')
        delete(aligned_data_file_rlb);
    end
    
    varNames = {'time','id','gas','signal'};

    ams_ts = this.get_from_binary( ams_fnames{1} ); % loading AMS TS with 1 sec frequency; assuming only one AMS TS
    
    Fts = this.sampl_res_const;
    Fact = this.deltaT_original;
    
    i_signal = 0;
    write_header = true;

    for i = 1:size(summary_table,1) % loop over the number of processed signals

        if ~summary_table(i,7) % only for reliable data
            continue;
        end
        
        if summary_table(i,6) ~= i_signal % only when we change the dataset
            i_signal = summary_table(i,6);
            snf_ts = this.get_from_binary(snf_fnames{i_signal,1}); % load Sniffer TS with 1 sec sampling frequency
        end

        first_ams = ( summary_table(i,2) - 1 ) * Fts + 1; % first index of the data range in ams_ts
        last_ams = ( summary_table(i,3) - 1 ) * Fts + 1; % last index of the data range in ams_ts
        first_snf = ( summary_table(i,4) - 1 ) * Fts + 1; % first index of the data range in snf_ts
        last_snf = ( summary_table(i,5) - 1 ) * Fts + 1; % last index of the data range in snf_ts

        % create resampling indeces
        l = 1;
        for j = first_ams:Fact:last_ams
            samples_ams(l,1) = uint32(j);
            l = l + 1;
        end
        l = 1;
        for j = first_snf:Fact:last_snf
            samples_snf(l,1) =  uint32(j);
            l = l + 1;
        end
    
        t = string( datetime( ams_ts( samples_ams,1 ),'ConvertFrom','datenum','Format','dd-MMM-yyy HH:mm:ss' ) );
        id = ams_ts( samples_ams,3 );
        gas = snf_ts(samples_snf,2:end);

        signal_no = zeros( numel(samples_snf),1 );
        signal_no(:,1) = i;

        if ( i > 1 )
            write_header = false;
        end

        aligned = table( t, id, gas, signal_no, 'VariableNames', varNames );

        writetable(aligned,aligned_data_file_rlb,'WriteMode','Append',...
                    'WriteVariableNames',write_header,'Delimiter',';',...
                    'QuoteStrings',true, 'FileType', 'text');

        clear samples_ams samples_snf signal_no
    end
   
end