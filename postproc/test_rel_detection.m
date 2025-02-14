clear;
close all

logs_table = read_logslist("logs_list.txt");

n_logs = height(logs_table);
data = cell(n_logs,1);
n_reliable = cell(n_logs,1);
for i = 1:n_logs
    data{i} = logs_table.files(i);
    n_reliable{i} = logs_table.nreliable(i);
end

% some notes
% 34: two-valued skew data
% 36: two-valued skew and very different data !!!
% 38: two-valued skew and very different data !!!
% 39: two-valued skew and very different data !!!
% 42: two-valued skew data, some reliable are negative
% 47: two-valued skew and very different data, one set with minuses !!!
% 48: two-valued skew and very different data !!!
% 51: two-valued skew and very different data
% 52: two-valued skew and very different data !!!
% 60: two-valued skew and very different data !!!
% 61: three-valued skew and very different data !!!
% 62: two-valued skew and very different data !!!

%% NOTE !!!

% This methods is not performing well for extreme cases where all data
% either reliable or not reliable at all. This is the limitation of the method,
% however, such cases need to be tasted and understood in order to be able to
% handle them.

%%
quality = 0.0;
range{1} = 2:64;
range{2} = 1:size(data,1);
range{3} = 38;
use_range = 2;

warning('off','all');
%warning('on','all');
%warning('query','all');

all_used_clusters = zeros( size(range{use_range},2), 3 );

for i = range{use_range}

    data_file = data{i};
    n_rlb = n_reliable{i};

    meth_opt_clst = 1; % use CalinskiHarabasz criterion

    [ all_res, opt_data, used_clusters ] = detect_reliable(data_file, meth_opt_clst);

    det_rlb = sum(all_res(:,2));
    accuracy = 1 - abs(det_rlb - n_rlb)/n_rlb;
    if accuracy < 0.0
        accuracy = 0.0;
    end
    quality = quality + accuracy;

    %disp([num2str(i), "clusters:", num2str(used_clusters), "detections:", num2str(det_rlb), "exp. detections:", num2str(n_rlb)]);
    disp( [ num2str(i), "accuracy, %:", num2str( accuracy ) ] );

    max_clusters = size(all_res,1); % max num clusters = data size
    all_used_clusters(i,1) = i;
    all_used_clusters(i,2) = used_clusters;
    all_used_clusters(i,3) = used_clusters/max_clusters;

    figure(1),
    clf;
    x = [min(opt_data(:,1)) max(opt_data(:,1))]; y = [n_rlb n_rlb];
    x2 = [used_clusters used_clusters]; y2 = [det_rlb det_rlb];

    subplot(2,1,1),
    yyaxis left; plot(opt_data(:,1),log(opt_data(:,2)), 'o-'); ylabel('log of optim. cluster measure')
    yyaxis right; plot(opt_data(:,1),opt_data(:,4), '*-'); ylabel('relative error change');
    xlabel('n.clusters / max.clusters')
    
    subplot(2,1,2),
    hold on, plot(opt_data(:,1),opt_data(:,3), 'o-'); plot(x,y); plot(x2,y2, 'p','MarkerFaceColor','red','MarkerSize',15);
    ylabel('num. detections'); xlabel('n.clusters / max.clusters')

end
disp([ "overall accuracy:", num2str(quality/size(data,1)) ])

%% Functions
function logslist = read_logslist(filename, dataLines)
    if nargin < 2
        dataLines = [1, Inf];
    end
    % Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 2);
    % Specify range and delimiter
    opts.DataLines = dataLines;
    opts.Delimiter = ";";
    % Specify column names and types
    opts.VariableNames = ["files", "nreliable"];
    opts.VariableTypes = ["string", "double"];
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    % Specify variable properties
    opts = setvaropts(opts, "files", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "files", "EmptyFieldRule", "auto");
    % Import the data
    logslist = readtable(filename, opts);
end