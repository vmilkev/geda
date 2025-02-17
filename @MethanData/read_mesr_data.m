function read_mesr_data( this )

this.MesrMap = containers.Map('KeyType','double','ValueType','any');
this.MesrNumMap = containers.Map('KeyType','double','ValueType','any');

data_sets = 0; % will be calculated based on detected gaps in data
data_sets_last = numel(this.geda_param.sniffile); % expected number of data sets

non_time_rec = 0;
short_time_rec = 0;

for i = 1:numel(this.geda_param.sniffile)

    [mesr, mesr_num] = this.importfile_csv2(this.geda_param.sniffile(i));
    
    tdiff = get_tdiff( mesr_num, this.signal_gap_limit );

    if ( isempty(tdiff) )
        tchange = mesr_num( 1:end, 1 ) - mesr_num( 1, 1 );
        if ( any( tchange < 0.0 ) )
            this.make_report("datint", "WARNING: Some Sniffer data records are not time-consecutive, in the file no.:", i);
            this.make_report("dat", "         This file will not be further processed!", []);
            this.make_report("dat", " ", []);
            continue;
        else
            data_sets = data_sets + 1;
            this.MesrMap( data_sets ) = mesr;
            this.MesrNumMap( data_sets ) = mesr_num;
            this.isMesrData = true;
        end
    else
        for j = 1:size(tdiff,1)
            tchange = mesr_num( tdiff(j,2):tdiff(j,3), 1 ) - mesr_num( tdiff(j,2), 1 );
            siglength = etime( datevec( mesr_num( tdiff(j,3), 1 ) ),datevec( mesr_num(tdiff(j,2), 1) ) )/( 60*60 ); % duration of sniffer records, [hours]
            if ( any( tchange < 0.0 ) || siglength < 1 )
                if ( any( tchange < 0.0 ) )
                    non_time_rec = non_time_rec + ( tdiff(j,3) - tdiff(j,2) );
                end
                if ( siglength < 1 )
                    short_time_rec = short_time_rec + ( tdiff(j,3) - tdiff(j,2) );
                end
            else
                data_sets = data_sets + 1;
                this.MesrMap( data_sets ) = mesr( tdiff(j,2):tdiff(j,3), : );
                this.MesrNumMap( data_sets ) = mesr_num( tdiff(j,2):tdiff(j,3), : );
            end
        end
        if ( tdiff(end,3) < size(mesr_num,1) )
            tchange = mesr_num( tdiff(end,3)+1:size(mesr_num,1), 1 ) - mesr_num( tdiff(end,3)+1, 1 );
            siglength = etime( datevec( mesr_num( size(mesr_num,1), 1 ) ),datevec( mesr_num(tdiff(end,3)+1, 1) ) )/( 60*60 ); % duration of sniffer records, [hours]
            if ( any( tchange < 0.0 ) || siglength < 1 )
                if ( any( tchange < 0.0 ) )
                    non_time_rec = non_time_rec + ( size(mesr_num,1) - tdiff(end,3)+1 );
                end
                if ( siglength < 1 )
                    short_time_rec = short_time_rec + ( size(mesr_num,1) - tdiff(end,3)+1 );
                end
            else
                data_sets = data_sets + 1;
                this.MesrMap( data_sets ) = mesr( tdiff(end,3)+1:size(mesr_num,1), : );
                this.MesrNumMap( data_sets ) = mesr_num( tdiff(end,3)+1:size(mesr_num,1), : );
            end
        end
        this.isMesrData = true;
    end
    if ( (data_sets - data_sets_last) > 0 )
        this.make_report("datint", strcat("WARNING! There are time gaps > ", num2str(this.signal_gap_limit), " sec present in the sniffer records. Data file no.:                  "), i);
        this.make_report("datint", "         The entire data will be divided into the number of datasets according to the discovered gaps: ", data_sets - data_sets_last );
        this.make_report("dat", " ", []);
        data_sets_last = data_sets;
    end
    if ( non_time_rec > 0 )
        this.make_report("datint", "WARNING: The records in some sub-datasets are not time-consecutive. The file no.:                           ", i);
        this.make_report("datint", "         In total, the following number of rows will not be further processed (due to time non continuity): ", non_time_rec );
        this.make_report("datint", "         The total number of records in the file:                                                           ", size(mesr_num,1) );
        this.make_report("dat", " ", []);
    end
    if ( short_time_rec > 0 )
        this.make_report("datint", "WARNING: The duration of some sub-dataset is less than 1 hour. The file no.:              ", i);
        this.make_report("datint", "         The following number of rows will not be further processed (due to time length): ", short_time_rec );
        this.make_report("datint", "         The total number of records in the file:                                         ", size(mesr_num,1) );
        this.make_report("dat", " ", []);
    end
end

    function diff = get_tdiff( data, gap ) % find gaps in a Sniffer signal
    t_diff = data(2:end,1) - data(1:end-1,1);
    lim_diff = gap/( 24*60*60 ); % allowed gap seconds to serial seconds
    i_d = find( t_diff >= lim_diff );
    if isempty(i_d)
        diff = [];
        return;
    end
    diff = zeros( size(i_d,1),3 );
    start = 1;
    for ii = 1:size(i_d,1)
        diff(ii,1) = t_diff( i_d(ii) );
        diff(ii,2) = start;
        diff(ii,3) = i_d(ii);
        start = i_d(ii) + 1;
    end
end

end