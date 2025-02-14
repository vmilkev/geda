function [res] = read_aligned(filename, dataLines)
    
    % reads resulting .txt file produced by GEDA into the table 'res'
    % example: [res] = read_aligned("aligned_init_sniffer_101.txt");

    if nargin < 2
        dataLines = [2, Inf];
    end

    opts = detectImportOptions(filename, 'FileType', 'text');

    % Specify range
    opts.DataLines = dataLines;
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    opts = setvaropts(opts, opts.VariableNames{1}, "InputFormat", "dd-MMM-yyyy HH:mm:ss");
    
    % Import the data
    res = readtable(filename, opts);

end