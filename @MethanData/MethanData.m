classdef MethanData < handle

    % PRIVATE INTERFACE
    properties (Access = private)
        % ................................................. LELY DATA:
        LelyMap % ......................................... Raw data, MAP format: KEYS <- robots, VALUES <-TABLE of data
        LelyNumMap % ...................................... Processed data, MAP format: KEYS <- robots, VALUES <-NUMERIC data
        Robots % .......................................... List of present robot names, array of NUMERIC format
        expected_robots % ................................. Expected names of robots, NUMERIC format
        isLelyData % ...................................... Flag indicating existance of Robot Data, LOGICAL format
        lely_param % ...................................... Data file parameters for LelyData, STRUCTURE format
        geda_param % ...................................... GEDA program parameters
        header_lely % ..................................... Names of columns of Lely NUMERIC data, array of STRINGS
        start_lely_stime % ................................ Very first serial day-time value in lely_ts data-set
        % ................................................. MESUREMENTS DATA:
        MesrMap % ......................................... Raw data, MAP format: KEYS <- file conseq. number, VALUES <-TABLE of data
        MesrNumMap % ...................................... Processed data, MAP format: KEYS <- file conseq. number, VALUES <-NUMERIC data
        isMesrData % ...................................... Flag indicating existance of Mesurements Data, LOGICAL format
        mesr_param % ...................................... Data file parameters for mesurements, STRUCTURE format
        co2_col % ......................................... Column in Mesr NUMERIC data pointing to CO2 record
        % ................................................. REPORT:
        text_fname % ...................................... File name with text report
        figs_fname % ...................................... File name with graphical report
        fig_resl % ........................................ Figures resolution if graphical report file
        f_hdl % ........................................... Report file handler
        sampl_res_const = 0 % ............................. Minimal resampling resolution of data, [sec]
        deltaT_original % ................................. Actual sampling frequency of the Sniffer data file
        signal_gap_limit = 120 % .......................... Minimal allowed sniffer signal gap, [sec]
    end
    properties (Access = public)
        device
    end

    % PUBLIK INTERFACE
    methods
        function this = MethanData() % .................... Class constructor
            % Set the defoult names
            this.expected_robots = [101 102 103];
            % We have not read data yet
            this.isLelyData = false;
            this.isMesrData = false;
            this.text_fname = [];
            this.figs_fname = [];
            this.fig_resl = '-r1200';
        end
        function delete( this ) % ......................... Class destructor
            fclose(this.f_hdl);
        end
        % ................................................. USER INTERFACE:
        set_robot_names( this, names ) % .................. Set expected robots' names
        set_lely_param( this, in_param ) % ................ Get necessary properties of LelyData file and set lely_param structure
        set_mesr_param( this, in_param ) % ................ Get necessary properties of MesrData file and set mesr_param structure
        [sampT,siglen,is_ssa,ssa1,ssa2] = set_geda_param(...
            this, in_param ) % Get necessary GEDA program parameters
        read_data( this ) % ............................... Reads .xlsx & .csv data files
        [data, robots] = get_lely_data( this, type ) % .... Public inteface to the class data members for Lely data
        data = get_mesr_data( this, type ) % .............. Public inteface to the class data members for Mesurements data
        [data, data_fnames] = get_lely_ts( this, robots ) % ............ Public inteface to the class data members for Lely TIME SERIES data
        [data, data2] = get_mesr_ts( this ) % ............. Public inteface to the class data members for Mesr TIME SERIES data
        r_name = get_mesr_robot( this ) % ................. Public inteface to the class data member this.mesr_param.robid
        set_report_param( this, t_name, f_name ) % ........ Set reporting parameters
        make_report( this, type, text, num_dat, ...
            map_dat, num_dat2, map_dat2 ) % ...... Publish report to the files
        res = adjust_skews( this, adj_table, ...
            lim_range, lim_time ) % ....... Finds right skews estimates and correct/adjust data output respectively.
        [ H,l ] = ssa2( this, X, L, Q ) % ................. Singular Spectrum Analysis implementation based on SVD, returns denoised data H
        c = cluster_data( this, data, n_clst, rep ) % ..... Clustering data using kmeans
        s = calc_cluster_statisstic( this, data, n_clst )
        set_sampling_freq( this, freq ) % ................. Assign value to the constaant sampl_res_const
        freq = get_sampling_freq( this ) % ................ Get value to the constant sampl_res_const
        freq = get_actsampling_freq( this ) % ............. Get value to the deltaT_original (actual sampling frequency)
        move_to_binary( this, data, fname ) % ............. Move array to a binary file
        data = get_from_binary( this, fname ) % ........... Extract saved data from a binary file
        remove_files(fname_list) % ........................ Delete binary files
        [ reliability_and_skew, optimization_data, n_clusters ] = detect_reliability( this, time_skew, use_matl_opt ) % Reliability evaluation based on estimated skews
        res = detect_reliability2( this, adj_table, lim_time ) % Second method (old) for reliability detection; works on smal datasets
        write_alignment_results( this, summary_table, ams_fnames, snf_fnames ) % Write alignment results into files
        fclean( this, files_list ) % remove *.bin files
        [res_stats, adj_tab, skew_res, test_stat1, test_stat2, ams_dt, snf_dt] = mf_detection(this, snf_data, ams_data, ams_arr, sig_len, use_ssa, report_msg) % Matched Filter signal detection
        [accuracy, corrected_rlb] = get_rlb_accuracy( this, rlb_estimates ) % Calculates the accuracy of data reliability estimation
    end

    % PRIVATE INTERFACE
    methods (Access = private)
        pstruct = read_fparam( this, filename ) % ......... Reads parameters' file
        read_lely_data(this) % ............................ Interface for importfile_xlsx()
        read_mesr_data(this) % ............................ Interface for importfile_csv()
        data = create_ts(this, in_data, header, delta_t) %  Transforms Lely numeric data to Time Series data
        [data, data2, d_range] = resample_ams(this, in_data) % ..... Transforms Lely numeric data to Time Series data
        data = scale( this, dat, col, flat, wndw ) % ...... Scalind data to the range [0,1]
        data = shrink( this, data, delta_t ) % ............ Function to reduce time interval between records
        [data1, data2, range, freq_mmed, freq_min, freq_max] = extend( this, d1, d2 ) % Function to extend time series to the range of 1 sec time interval between records
        [data1, data2, range, freq_mmed, freq_min, freq_max] = resample_snf( this, d1, d2 )
        [mesr, mesr_num] = importfile_csv( ... % .......... Reads .csv file
            this, ...
            filename, ...
            dataTypes, ...
            dataLines)
        [mesr, mesr_num] = importfile_csv2(this, filename) % Reads .csv file, only Sniffer data
        [map, num_map, robots] = importfile_xlsx( ... % ... Reads .xlsx file
            this, ...
            workbookFile, ...
            dataTypes, ...
            sheetName, ...
            dataLines)
        [map, num_map, robots] = importfile_xlsx2(this) % ... Reads .xlsx file, only AMS data
    end

end

