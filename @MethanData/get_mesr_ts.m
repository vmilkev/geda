function [data, data_noscale] = get_mesr_ts( this, delta_t )

% PUBLIK INTERFACE TO THE CLASS MEMBER VARIABLES
% RETURNS A MESUREMENTS TIME SERIES DATA

if ~this.isMesrData
    this.read_data();
end

if isempty(this.MesrMap)
    disp("get_mesr_ts(): The Sniffer data is empty!");
    return;
end

if delta_t < 1
    disp("get_mesr_ts(): The sampling time interval (deltaT) should be higher than 1 sec.");
    return;
end

data = cell( length(this.MesrMap),1 );
data_noscale = cell( length(this.MesrMap),1 );

this.co2_col = 3;%find(this.mesr_param.head == "co2");
tcols = 1;%find(this.mesr_param.head == "time");
%this.co2_col = this.co2_col - ( numel(tcols) - 1);

if numel(tcols) == 0
    disp("get_mesr_ts(): There is no TIME record in the data.");
    return;
end

if isempty(this.MesrMap)
    disp("get_mesr_ts(): this.MesrMap is empty.");
    return;
end

if isempty(this.MesrNumMap)
    disp("get_mesr_ts(): this.MesrNumMap is empty.");
    return;
end

pct = 50; % use a specific percentile

for i = 1:length(this.MesrMap)

    d = this.scale( this.MesrNumMap(i), this.co2_col, true, pct );
    d2 = this.MesrNumMap(i); % no scaling
    
    % extend tata to 1 sec sampling interval
    [d,d2] = this.extend( d, d2 );

    if delta_t == 1
        data{i,1} = d;
        data_noscale{i,1} = d2;
        continue;
    end
    if nargin < 2
        % default shrinkage to 5 sec interval
        data{i,1} = this.shrink( d );
        data_noscale{i,1} = this.shrink( d2 );
    else
        data{i,1} = this.shrink( d, delta_t );
        data_noscale{i,1} = this.shrink( d2, delta_t );
    end
end

end