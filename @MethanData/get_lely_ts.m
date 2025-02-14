function [data, data_fnames] = get_lely_ts( this, robots )

% PUBLIK INTERFACE TO THE CLASS MEMBER VARIABLES
% RETURNS AMS TIME SERIES DATA

if ~this.isLelyData
    this.read_data();
end

if isempty(this.Robots)
    disp("get_lely_ts(): The ROBOT data is empty.");
    return;
end

data = cell( numel(robots),1 );
data_stat = zeros(numel( robots ),4);

for i = 1:numel( robots )

    tic;

    [data{i,1}, data2, data_range] = this.resample_ams( this.LelyNumMap(robots(i)) ); % resample to delta_t interval

    fname = strcat("ams_",num2str(i),".bin");
    this.move_to_binary(data2,fname);
    clear data2
    data_fnames{i,1} = fname;

    elps_time = toc;

    data_stat(i,1) = i; % data set no.
    data_stat(i,2) = floor(data_range*24); % data length, hours
    data_stat(i,3) = this.get_sampling_freq(); % resampling freq, sec
    data_stat(i,4) = elps_time; % elapsed time

end

this.make_report("stats", "AMS data summary:", data_stat);

end