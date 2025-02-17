function [lely_map, lely_num_map, robots] = importfile_xlsx2(this)

% Set up the Import Options and import the data

lines = this.geda_param.amsrange;
workbookFile = this.geda_param.amsfile;

opts = detectImportOptions(workbookFile);
opts.DataLines = lines;

opts.SelectedVariableNames = {char(this.geda_param.id),...
                              char(this.geda_param.time(1,1)),...
                              char(this.geda_param.time(2,1))...
                              };
LelyData = readtable(workbookFile, opts);

opts.SelectedVariableNames = {char(this.geda_param.device)};
device = table2array( readtable(workbookFile, opts) );

% Create LELYMAP object & LELY_NUM_MAP object
%
% Create a map object where keys are robot names and values are arrays of
% matrix data for a specific robot

% Define empty Map
lely_map = containers.Map('KeyType','double','ValueType','any');

% Define empty Map
lely_num_map = containers.Map('KeyType','double','ValueType','any');

% Array of unique robot names present in the data
robots = this.get_mesr_robot(); % do not extract all robots, extract only the one which corresponds to sniffer data file

for i = 1:numel(robots)
    [rows,~] = find( device(:,1) == robots(i) );
    lely_map( robots(i) ) = LelyData(rows,:);
    
    % Transform lely_map to numerical data map
    arr = zeros( numel(rows), 3 );
    arr(:,1) = LelyData{rows,1};

    if isa(LelyData{1,2},'datetime')
        arr(:,2:3) = datenum( LelyData{rows,2:3} );
    else
        dt2 = LelyData{rows,2};
        dt3 = LelyData{rows,3};

        dt2 = strrep(dt2,"'",'');
        dt3 = strrep(dt3,"'",'');

        dt2 = datetime(dt2);
        dt3 = datetime(dt3);

        %fmt = 'yyyy-mm-dd HH:MM:SS';
    
        arr(:,2) = datenum( dt2 );%, fmt );
        arr(:,3) = datenum( dt3 );%, fmt );
    end
    inan_2 = isnan( arr(:,2) );
    inan_3 = isnan( arr(:,3) );
    inan = inan_2 | inan_3;

    n_nan = numel( find(inan == 1));
    if ( n_nan > 1 )
        this.make_report( "dat", "WARNING: NAN records were found in the AMS data file!                               ", [] );
        this.make_report( "dat", "         The following number of records will be excluded from the processing data: ", n_nan );
        this.make_report( "dat", "         The number of records will be processed:                                   ", numel(inan) - n_nan );
        this.make_report( "dat", " ", []);
    end

    lely_num_map( robots(i) ) = arr(~inan,:);

    for i2 = 1:numel(inan)
        if ( inan(i2) == false )
            this.start_lely_stime = arr(i2,2);
            break;
        end
    end

end

end