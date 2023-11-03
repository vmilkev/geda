function [lely_map, lely_num_map, robots] = importfile_xlsx(this, workbookFile, dataTypes, sheetName, dataLines)

%%
% IMPORT DATA FROM A SPREADSHEET
%
%   MAP is a map object where Keys are robot names and values are a
%           data records for a specific robot;
%
%   ROBOTS is an array containing robot names presented in a data file;
%
%   NUM_MAP is a map object where Keys are robot names and values are
%           a numerical data records for a specific robot;
%
%  [MAP, NUM_MAP, ROBOTS] = IMPORTFILE_XLSX(FILE, TYPES)
%             
%             Reads from the specified FILE for all rows, in a sheet 1.
%             Specify TYPES as a string array corresponding to data
%             columns, example: ["double", "string", "datetime"].
%
%  [MAP, NUM_MAP, ROBOTS] = IMPORTFILE_XLSX(FILE, TYPES, SHEET)
%             
%             Reads from the specified FILE for all rows.
%             Specify TYPES as a string array corresponding to data
%             columns, example: ["double", "string", "datetime"].
%             Specify a SHEET name as a string, example: 'Ark1'.
%
%  [MAP, NUM_MAP, ROBOTS] = IMPORTFILE_XLSX(FILE, TYPES, SHEET, LINES)
%             
%             Reads from the specified FILE for all rows.
%             Specify TYPES as a string array corresponding to data
%             columns, example: ["double", "string", "datetime"].
%             Specify a SHEET name as a string, example: 'Ark1'.
%             Specify a rows range by LINES as an integer array, example: [2, 10].
%
%
%  EXAMPLE:
%
%  [map, num_map, robots] = importfile_xlsx(...
%             "/data/LelyData_21834.xlsx", ...
%             ["double", "string", "double", "double", "double", "datetime", "datetime"], ...
%             "Ark1", ...
%             [2, Inf]);
%

%% Input handling

if nargin < 2
    disp("There is not enough parameters. Exit.");
    return;
end

% If dataLines & sheetName are not specified, define defaults
if nargin < 3
    sheetName = 1;
    dataLines = [2, Inf];
end

% If row start and end points are not specified, define defaults
if nargin < 4
    dataLines = [2, Inf];
end

%% Determine file type

is_xlsx = false;
is_csv = false;

if contains( workbookFile, '.xlsx' )
    is_xlsx = true;
elseif contains( workbookFile, '.csv' )
    is_csv = true;
else
    disp(["importfile_xlsx(): Wrong type of data file ", workbookFile, " Expected types: xlsx or csv"]);
    return;
end

%% Set up the Import Options and import the data

if is_xlsx
    
    opts = spreadsheetImportOptions("NumVariables", numel(dataTypes));

    % Specify sheet and range
    opts.Sheet = sheetName;
    opts.DataRange = dataLines;

    % Specify column types
    opts.VariableTypes = dataTypes; %["double", "string", "double", "double", "double", "datetime", "datetime"];

    % Import the data
    LelyData = readtable(workbookFile, opts, "UseExcel", false);
    
else
    
    opts = detectImportOptions(workbookFile);
    opts.DataLines = dataLines;

    if strcmp(this.lely_param.ams, "delaval")
%         opts = detectImportOptions(workbookFile);
%         opts.DataLines = dataLines;
        opts.SelectedVariableNames = {'DeviceName'};
        delval_rbt_name = readtable(workbookFile, opts); %delval_rbt_name should be "DeviceName";
    end

    opts = delimitedTextImportOptions("NumVariables", numel(dataTypes));

    % Specify range and delimiter
    %opts.DataLines = dataLines;

%     if strcmp(this.lely_param.ams, "lely")
%         opts.Delimiter = ";";
%     else
%         opts.Delimiter = ",";
%     end

    opts.VariableTypes = dataTypes;

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    LelyData = readtable(workbookFile, opts);
    
end
%% Create LELYMAP object & LELY_NUM_MAP object
%
% Create a map object where keys are robot names and values are arrays of
% matrix data for a specific robot

% Define empty Map
lely_map = containers.Map('KeyType','double','ValueType','any');

% Define empty Map
lely_num_map = containers.Map('KeyType','double','ValueType','any');

% list of expected robots
%expected_robots = [101 102 103];

% Find possible cols. where the robot names could be
numCols = find(dataTypes == "double");
numCols2 = find(dataTypes == "double" | dataTypes == "datetime");
numCols3 = find(dataTypes == "datetime");

% get the header of numeric data
this.header_lely = this.lely_param.head(numCols2);

% Find the column where robot names are recorded
for j2 = 1:numel(this.expected_robots)
    for i2 = 1:numel(numCols)
        [~,l2] = find(LelyData{:,numCols(i2)} == this.expected_robots(j2));
        if ( ~isempty(l2) )
            break;
        end
    end
end

% Array of unique robot names present in the data
robots = table2array( unique( LelyData(:,numCols(i2)) ) );

for i = 1:numel(robots)
    [rows,~] = find(LelyData{:,numCols(i2)} == robots(i));
    lely_map( robots(i) ) = LelyData(rows,:);
    
    % Transform lely_map to numerical data map
    arr = zeros( numel(rows), numel(numCols2) );
    arr(:,1:numel(numCols)) = LelyData{rows,numCols};
    arr(:,numel(numCols)+1:numel(numCols)+numel(numCols3)) = datenum( LelyData{rows,numCols3} );   
    lely_num_map( robots(i) ) = arr;
end

end