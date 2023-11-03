function [mesr, mesr_num] = importfile_csv(this, filename, dataTypes, dataLines)
%
% IMPORT DATA FROM CSV FILE
%
%  MESR is a map with table data
% 
%  MESR_NUM is a map with numeric data
% 
%  [MESR, MESR_NUM] = IMPORTFILE_CSV(FILE)
%             
%             Reads from the specified FILE for all rows.
%             As defoult, TYPES = ["categorical", "datetime", "double",
%             "double"]; and LINES = [2, Inf].
%
%  [MESR, MESR_NUM] = IMPORTFILE_CSV(FILE, TYPES)
%             
%             Reads from the specified FILE for all rows.
%             Specify TYPES as a string array corresponding to data
%             columns, example: ["double", "string", "datetime"].
%
%  [MESR, MESR_NUM] = IMPORTFILE_CSV(FILE, TYPES, LINES)
%             
%             Reads from the specified FILE for all rows.
%             Specify TYPES as a string array corresponding to data
%             columns, example: ["double", "string", "datetime"].
%             Specify a rows range by LINES as an integer array, example: [2, 10].
%
%  EXAMPLE:
%
%  [mesr, mesr_num] = importfile_csv(...
%             "data/21834_001_16102018_1.csv", ...
%             ["categorical", "datetime", "double", "double", "string", "string"], ...
%             [2, Inf]);
%
%% Input handling

if nargin < 1
    disp("importfile_csv(): There is not enough parameters. Exit.");
    return;
end

% If dataTypes & dataLines are not specified, define defaults
if nargin < 2
    dataTypes = ["categorical", "datetime", "double", "double"];
    dataLines = [2, Inf];
end

% If dataLines is not specified, define defaults
if nargin < 3
    dataLines = [2, Inf];
end

%%
% Change comma to dot in numbers
Data = fileread(filename);
Data = strrep(Data, ',', '.');
tFile = strcat(filename,'2');
FID = fopen( tFile, 'w');

% command = 'pwd';
% [status,cmdout] = system(command)

fwrite(FID, Data, 'char');
fclose(FID);
filename = tFile;

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", numel(dataTypes));

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

opts.VariableTypes = dataTypes; %["categorical", "datetime", "double", "double", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
mesr = readtable(filename, opts);
sz2 = size(mesr);

if isempty(mesr)
    disp("importfile_csv(): readtable(filename, opts) gives empty mesr variable.");
    return;
end

%%
% Find possible cols. where the robot names could be
numCols = find(dataTypes == "categorical" | dataTypes == "datetime");
numCols2 = find(dataTypes == "double");

% Modify day/time of robot record robRec_0 and write it to the new array robRec
mesr_num = zeros( sz2(1,1), numel(numCols2)+1 );
fmt = 'dd-mm-yyyy HH:MM:SS';

d = cellstr(mesr{:,numCols(1)});
t = cellstr(mesr{:,numCols(2)});
dt = strcat( d, {' '}, t );

this.make_report("dat", "importfile_csv(): converting day-time using format dd-mm-yyyy HH:MM:SS", []);
this.make_report("dat", "importfile_csv(): in the case of error, check a day-time format in a data file!", []);

mesr_num( :,1 ) = datenum( dt,fmt );
mesr_num( :,2:numel(numCols2)+1 ) = mesr{:,numCols2};

end