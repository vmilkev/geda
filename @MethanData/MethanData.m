classdef MethanData < handle

    % PRIVATE INTERFACE
    properties (Access = private)
        % ................................................. LELY DATA:
        LelyMap % ......................................... Raw data, MAP format: KEYS <- robots, VALUES <-TABLE of data
        LelyNumMap % ...................................... Processed data, MAP format: KEYS <- robots, VALUES <-NUMERIC data
        Robots % .......................................... List of present robot names, array of NUMERIC format
        isLelyData % ...................................... Flag indicating existance of Robot Data, LOGICAL format
        geda_param % ...................................... GEDA program parameters
        header_lely % ..................................... Names of columns of Lely NUMERIC data, array of STRINGS
        start_lely_stime % ................................ Very first serial day-time value in lely_ts data-set
        % ................................................. MESUREMENTS DATA:
        MesrMap % ......................................... Raw data, MAP format: KEYS <- file conseq. number, VALUES <-TABLE of data
        MesrNumMap % ...................................... Processed data, MAP format: KEYS <- file conseq. number, VALUES <-NUMERIC data
        isMesrData % ...................................... Flag indicating existance of Mesurements Data, LOGICAL format
        co2_col % ......................................... Column in Mesr NUMERIC data pointing to CO2 record
        % ................................................. REPORT:
        text_fname % ...................................... File name with text report
        f_hdl % ........................................... Report file handler
        % .................................................
        sampl_res_const = 0 % ............................. Minimal resampling resolution of data, [sec]
        downsampl_res = 0 % ............................... Resampling resolution of data used for study the effect of downsampling, [sec]
        deltaT_original % ................................. Actual sampling frequency of the Sniffer data file
        signal_gap_limit = 120 % .......................... Minimal allowed sniffer signal gap, [sec]
        % .................................................
        outpath % ......................................... Output directory for resulting files
        % ................................................. Outliers filtering:
        bkg_outl_const = -0.1 % ........................... background outlier upper bound
        obs_outl_const = 1.5 % ............................ observations outlier in statistic: data./(3*std(data)) > obs_outl_const => is outlier (if expression is TRUE);
        wlen_const = 10 % ................................. SSA window for denoising, sec
        explvar_const = 90 % .............................. amount of variance a denoised data should cover/explain
    end
    properties (Access = public)
        device % robot id, assumed numeric type (double) !!!
        sniffer = "none" % sniffer id, assume striung
        farm = "none" % farm id, assume string
        is_ssa = true;
        ssa_wnd = 200;
        ssa_eigval = 5;
        date_now = ""; % day-time without sopaces and colomns in between
        % .................................................
        trait_type = 1 % .................................. Determines which type of trait to output: 1 (default) => time averagined value; 2 => averaged over num of observations
    end

    % PUBLIK INTERFACE
    methods
        function this = MethanData() % .................... Class constructor
            this.isLelyData = false;
            this.isMesrData = false;
            this.text_fname = [];
            this.date_now = regexprep(regexprep(string(datetime()), ' +', '_'), ':+', '-');
        end
        function delete( this ) % ......................... Class destructor
        end
        function on_exit( this ) % ......................... Instead of Class destructor, since we have problem in parfor loops on class copies
            fclose(this.f_hdl); % closing the log file
        end
        % ................................................. USER INTERFACE:
        set_robot_names( this, names ) % .................. Set expected robots' names
        [sampT,siglen] = set_geda_param(this, in_param ) %  Get necessary GEDA program parameters
        read_data( this ) % ............................... Reads .xlsx & .csv data files
        [data, robots] = get_lely_data( this, type ) % .... Public inteface to the class data members for Lely data
        data = get_mesr_data( this, type ) % .............. Public inteface to the class data members for Mesurements data
        [data, data_fnames] = get_lely_ts( this, robots ) % Public inteface to the class data members for Lely TIME SERIES data
        [data, data2] = get_mesr_ts( this ) % ............. Public inteface to the class data members for Mesr TIME SERIES data
        r_name = get_mesr_robot( this ) % ................. Public inteface to the class data member this.mesr_param.robid
        set_report_param( this, t_name ) % ................ Set reporting parameters
        make_report( this, type, text, num_dat, map_dat ) % Publish report to the files
        [ H,l,nl,F ] = ssa( this, X, L, Q, E ) % .......... Singular Spectrum Analysis implementation based on Covariance matrix and PCA, returns denoised data H
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
        fclean( this, files_list ) % ...................... remove *.bin files
        [res_stats, adj_tab, skew_res, test_stat1, test_stat2, ams_dt, snf_dt] = mf_detection(this, snf_data, ams_data, ams_arr, sig_len, use_ssa, report_msg) % Matched Filter signal detection
        [accuracy, corrected_rlb] = get_rlb_accuracy( this, rlb_estimates ) % Calculates the accuracy of data reliability estimation
        res = apply_ssa_ams(this, arr, max_length, ssa_window, ssa_eigenvalue)
        write_traits( this, summary_table, ams_fnames, snf_fnames, used_sig_length )
    end

    % PRIVATE INTERFACE
    methods (Access = private)
        pstruct = read_fparam( this, filename ) % ......... Reads parameters' file
        read_lely_data(this) % ............................ Interface for importfile_xlsx()
        read_mesr_data(this) % ............................ Interface for importfile_csv()
        [data, data2, d_range] = resample_ams(this, in_data) % ..... Transforms Lely numeric data to Time Series data
        data = scale( this, dat, col, flat, wndw ) % ...... Scalind data to the range [0,1]
        [data1, data2, range, freq_mmed, freq_min, freq_max] = resample_snf( this, d1, d2 )
        [mesr, mesr_num] = importfile_csv2(this, filename) % Reads .csv file, only Sniffer data
        [map, num_map, robots] = importfile_xlsx2(this) % .. Reads .xlsx file, only AMS data
    end

end

