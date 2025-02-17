function [mesr, mesr_num] = importfile_csv2(this, filename)

%this.geda_param.range
% Set up the Import Options and import the data
opts = detectImportOptions(filename);

for i = 1:numel(opts.VariableTypes)
    if ( strcmp(opts.VariableTypes{1,i},'double') )
        opts = setvartype(opts,{opts.VariableNames{1,i}},'char');
    end
end

% Specify range
opts.DataLines = this.geda_param.snifrange;

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
mesr = readtable(filename, opts);

if isempty(mesr)
    this.make_report( "dat", "ERROR in importfile_csv2(): readtable(filename, opts) gives empty mesr variable! ", [] );
    return;
end

sz2 = size(mesr);

% Modify day/time of robot record robRec_0 and write it to the new array robRec
mesr_num = zeros( sz2(1,1), sz2(1,2)-1 );

d = cellstr(mesr{:,1});
t = cellstr(mesr{:,2});
dt = strcat( d, {' '}, t );

dt = strrep(dt,'/', '-');
dt = strrep(dt,'.', ':');

% select apropriate format
fmt{1} = 'dd-mm-yyyy HH:MM:SS'; % default format accepted from the start
fmt{2} = 'yyyy-mm-dd HH:MM:SS';
%fmt{3} = 'yyyy-MM-dd HH:mm:SS';
%fmt{4} = 'dd/MM/yyyy HH:mm:SS';
%fmt{4} = 'dd.MM.yyyy HH:mm:SS';
%fmt{5} = 'yyyy/MM/dd HH:mm:SS';
%fmt{6} = 'yyyy.MM.dd HH:mm:SS';

num_formats = numel(fmt);

time_diff = zeros(num_formats,1);

for i_fmt = 1:num_formats
    time_diff(i_fmt,1) = abs(this.start_lely_stime - datenum( dt{1,1},fmt{i_fmt} ))*24*60*60; % diff is seconds
end
[~,min_fmt] = min( time_diff );

%datetime.setDefaultFormats('default',fmt{min_fmt}); % setting the default format

this.make_report("dat",  strcat( "Converting Sniffers day-time records using the following format: ",convertCharsToStrings( fmt{min_fmt} ) ), []);
this.make_report("dat", "NOTE: Check if a day-time format in a Sniffer data file corresponds to the one used in the program.", []);
this.make_report("dat", " ", []);

mesr_num( :,1 ) = datenum( dt,fmt{min_fmt} );
%mesr_num( :,1 ) = datenum( dt ); % here we are using the default format

for i_col = 3:sz2(1,2)
    i_val = mesr{1,i_col};
    if isa(i_val, 'cell')
        val = mesr{:,i_col};
        val = strrep(val,',', '.');
        mesr_num( :,i_col-1 ) = str2double(val);
    elseif isa(i_val, 'datetime')
        this.make_report("dat", "WARNING: In the sniffer data file a column consists of datetime data, will be omitted. Col no.: ", i_col);
        this.make_report("dat", " ", []);
        continue;
    elseif isa(i_val, 'double')
        mesr_num( :,i_col-1 ) = mesr{:,i_col};
    else
        this.make_report("dat", "WARNING: undefined data type found in a column of the sniffer data file. Col no.: ", i_col);
        this.make_report("dat", " ", []);
        return;
    end
end

inan_2 = isnan( mesr_num(:,2) );
inan_3 = isnan( mesr_num(:,3) );
inan = inan_2 | inan_3;

n_nan = numel( find(inan == 1));
if ( n_nan > 1 )
    this.make_report( "dat", "WORNING: NAN records were found in the sniffer data file.                           ", [] );
    this.make_report( "dat", "         The following number of records will be excluded from the processing data: ", n_nan );
    this.make_report( "dat", "         The following number of records will be processed:                         ", numel(inan) - n_nan );
    this.make_report( "dat", " ", []);
end

mesr_num = mesr_num(~inan,:);
mesr = mesr(~inan,:);

end