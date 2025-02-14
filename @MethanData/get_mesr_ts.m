function [data,data_noscale] = get_mesr_ts( this )

% PUBLIK INTERFACE TO THE CLASS MEMBER VARIABLES
% RETURNS A MESUREMENTS TIME SERIES DATA

if ~this.isMesrData
    this.read_data();
end

if isempty(this.MesrMap)
    warning("get_mesr_ts(): The Sniffer data is empty!");
    return;
end

data = cell( length(this.MesrMap),1 );
data_noscale = cell( length(this.MesrMap),1 );

this.co2_col = 3;%find(this.mesr_param.head == "co2");
tcols = 1;%find(this.mesr_param.head == "time");

if numel(tcols) == 0
    this.make_report("dat", "WARNING in get_mesr_ts(): There is no TIME records in the data!", []);
    return;
end

if isempty(this.MesrMap)
    this.make_report("dat", "WARNING in get_mesr_ts(): this.MesrMap is empty!", []);
    return;
end

if isempty(this.MesrNumMap)
    this.make_report("dat", "WARNING in get_mesr_ts(): this.MesrNumMap is empty!", []);
    return;
end

pct = 50; % use a specific percentile

data_stats = zeros(length(this.MesrMap),6);

for i = 1:length(this.MesrMap)
    
    tic;

    d = this.scale( this.MesrNumMap(i), this.co2_col, true, pct ); % rescaling data
    d2 = this.MesrNumMap(i); % no rescaling
    
    [d,d2,d_range,freq_median,freq_min,freq_max] = this.resample_snf( d, d2 ); % resampling tata to the native (median) sampling rate
    
    fname = strcat("snf_",num2str(i),".bin");
    this.move_to_binary(d2,fname);
    clear d2

    data{i,1} = d;
    data_noscale{i,1} = fname;

    elps_time = toc;

    data_stats(i,1) = i; % Sniffer data set no.
    data_stats(i,2) = floor(d_range); % Sniffer data range (between start & end of records), whole hours
    data_stats(i,3) = int32(freq_median); % MEDIAN of sampling frequency, sec
    data_stats(i,4) = int32(freq_min); % MIN of sampling frequency, sec
    data_stats(i,5) = int32(freq_max); % MAX of sampling frequency, sec
    data_stats(i,6) = elps_time; % Elapsed time for processing

end

this.make_report("stats", "Sniffer data summary:", data_stats);

end