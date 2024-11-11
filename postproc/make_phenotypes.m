function make_phenotypes( fname, w_length, filter_ids_fname )

% Takes as an input the GEDA produced file
% (such as aligned_init_sniffer_101.txt)
% and do the following (for each gas type present):
% 1. calculates the backgroung values;
% 2. correct data for the background
% 3. do SSA-based denoising !!! Not in current version
% 4. calculates emission phenotype.
%
% fname - input file name
% w_length - max duration in hours used for background calculations

if nargin < 1
    disp("At least a GEDA resulting file with aligned data is expected as an input argument!");
    return;
end

if nargin < 3
    if ~exist('w_length','var')
      w_length = 6;
    end
    if ~exist('fname','var')
        disp("At least a GEDA resulting file with aligned data is expected as an input argument!");
        return;
    end
end

data = read_aligned(fname); % read file

signal_id_list = [];

if exist('filter_ids_fname','var')
    signal_id_list = read_signals_ids(filter_ids_fname);

    if ~isempty(signal_id_list)
        data = filter_signals(data, signal_id_list);
        %writetable(data,"filtered_data.txt",'WriteMode','overwrite');
        %return;
    end
end

t_dif2 = data.time(2:end, 1)-data.time(1:end-1,1);
t_mask = t_dif2 < 0;
if sum(t_mask)/numel(t_mask) < 0.0001
    data(t_mask,:) = [];
    t_dif2 = data.time(2:end, 1)-data.time(1:end-1,1);
    t_mask = t_dif2 < 0;
    while sum(t_mask) > 0
        data(t_mask,:) = [];
        t_dif2 = data.time(2:end, 1)-data.time(1:end-1,1);
        t_mask = t_dif2 < 0;
        if sum(t_mask)/numel(t_mask) > 0.0001
            disp("The time in the data set is not consecutive. The data needs to be revised.");
            return;
        end
    end
else
    disp("The time in the data set is not consecutive. The data needs to be revised.");
    return;
end

time_col = 1;
id_col = 2;
outlier_low_limit = -0.1;

% determine the number of gas types
res_dim = size(data);
n_gases = res_dim(1,2) - 3; % 3 because: time, id, and signal cols
max_data_rows = res_dim(1,1);

range0 = 1:max_data_rows; % use the entire range

%range0 = 1:floor(max_data_rows/4); % for testing/debugging use a part of the data
%range0 = floor(max_data_rows/4)+1:floor(max_data_rows/4)*2;
%range0 = floor(max_data_rows/4)*2+1:floor(max_data_rows/4)*3;
%range0 = floor(max_data_rows/4)*3+1:floor(max_data_rows/4)*4;

id = table2array(data(range0,id_col)); % convert tables to arrays
t = convertTo(table2array(data(range0,time_col)), 'posixtime'); % convert tables to arrays

for igas = (id_col+1):res_dim(1,2)-1 % loop over gas types
    
    gas = table2array(data(range0,igas)); % get gas data, convert tables to arrays

    if std(gas) == 0
        continue;
    end

    if sum(isnan(gas)) > 0
        continue;
    end
    
    mask = find(gas >= outlier_low_limit); % mask outliers
        
    d = [t(mask) id(mask) gas(mask)]; % all data without filtered outliers
    
    mask = d(:,id_col) == 0; % mask the background dataset
    
    d0 = d(mask,:); % background dataset
    d1 = d(~mask,:); % not-background dataset

    clear d gas
    
    day_in_sec = 24*60*60; % duration of one day, seconds
    n_days = ceil( ( d0(end, 1) - d0(1,1) ) / day_in_sec ); % number of days in the dataset, days
    wnd_len = min(w_length*60*60, d0(end, 1) - d0(1,1)); % duration of moving window, seconds
    
    time_pos  = find_pos_by_value( d0(:,1), wnd_len );
    
    n_wnd = size(time_pos,1)+1; % num of windows for the entire record
    day_periods = ceil(n_wnd/n_days); % num of windows per day
       
    bkg = get_background( d0, time_pos, n_wnd ); % calculate background concentration for each window
    pos_arr = find_pos_by_value2( bkg(:,1), d1(:,1) );
    
    % plot_overview();
    % plot_bkg();
    % plot_bkg_obs();
    % plot_bkg_uncorr_obs();

    for i = 1:size(pos_arr,1)-1 % correct for backkground
        d1(pos_arr(i):pos_arr(i+1)-1,3) = d1(pos_arr(i):pos_arr(i+1)-1,3) - bkg(i,2);
    end
    d1(pos_arr(end):end,3) = d1(pos_arr(end):end,3) - bkg(end,2);
    mask = d1(:,3) > 0; % find negative outliers
    d1 = d1(mask,:); % exclude negative outliers

    clear d0 bkg pos_arr time_pos mask

    unique_ids = unique( d1(:,2) ); % get list of ids
    %traits = cell( size(unique_ids,1), size(unique_ids,1) ); % each cell will consist results for a specific id
    traits = zeros(ceil(size(unique_ids,1) * 10), 6); % each cell will consist results for a specific id
    j = 1;
    for i = 1:size(unique_ids,1) % loop over specific ids
        mask = d1(:,2) == unique_ids(i); % get all records for the specific id
        d_id = d1(mask,:); % make copy for convenience
        t_gaps = d_id(2:end,1) - d_id(1:end-1,1); % get records difference in order to find big gaps
        md = get_mode(t_gaps); % find sampling interval/frequency
        gaps = find( t_gaps > md * 10 ); % find time gaps between the measurements of the specific id
        trt = zeros(1,6);
        if ~isempty(gaps) % if there are more than one set of measurements
            r1 = 1;
            for h = 1:size(gaps,1) % loop over distant records of a specific id
                i_range = r1:gaps(h);
                if h == size(gaps,1)
                    i_range = (gaps(h)+1):size(d_id,1);
                end
                trt(1,1:2) = d_id(i_range(1), 1:2); % time and id
                trt(1,3) = get_trait( d_id(i_range, :), 0.1, "integral" );
                trt(1,4) = get_trait( d_id(i_range, :), 0.1, "peaks" );
                trt(1,5) = get_trait( d_id(i_range, :), 0.1, "peaksintegral" );
                trt(1,6) = get_trait( d_id(i_range, :), 0.1, "average" );

                traits(j, :) = trt;

                %traits = [ traits; d_id(i_range, :) ];
                %traits{i,h} = d_id(i_range, :);
                
                r1 = gaps(h) + 1;
                j = j + 1;
            end
        else % there is only one set of records
            trt(1,1:2) = d_id(1, 1:2); % time and id
            trt(1,3) = get_trait( d_id, 0.1, "integral" );
            trt(1,4) = get_trait( d_id, 0.1, "peaks" );
            trt(1,5) = get_trait( d_id, 0.1, "peaksintegral" );
            trt(1,6) = get_trait( d_id, 0.1, "average" );

            traits(j, :) = trt;

            j = j + 1;
            %traits{i,1} = d_id;
            %traits = [ traits; d_id ];
        end
    end

    mask = traits(:,3) > 0;
    traits = traits(mask, :);

    varNames = {'time','id','trait_1','trait_2','trait_3','trait_4'};
    results = table( datetime(traits(:,1), 'convertfrom', 'posixtime'), traits(:,2), traits(:,3), traits(:,4), traits(:,5), traits(:,6),'VariableNames',varNames );
    
    k = strfind(fname, '/');
    if ~isempty(k)
        fname2 = extractAfter(fname,k(end)+1);
    else
        fname2 = fname;
    end

    gas_trait_file = strcat("traits_gas_", num2str(igas-2), "_", fname2);
    writetable(results,gas_trait_file,'WriteMode','overwrite',...
            'WriteVariableNames',true,'Delimiter',';','QuoteStrings',true);


end % end of for_loop ove gas types

% NESTED FUNCTIIONS

    function plot_overview()
        figure(1);
        subplot(2,1,1), plot(d0(:,1)./60, d0(:,3), 'bo');
        subplot(2,1,2), plot(d1(:,1)./60, d1(:,3), 'ro');
    end

    function plot_bkg()
        figure(2);
        t0 = d0(1,1);
        plot(bkg(:,1)./60 -t0./60, bkg(:,2), '-o');
    end

    function plot_bkg_obs()
        figure(3)
        r2 = 0;
        t0 = d0(1,1);
        for i = 1:n_wnd-1
            r1 = r2 + 1;
            r2 = time_pos(i);
            subplot(n_days, day_periods, i);
            plot(d0(r1:r2,1)./60-t0./60, d0(r1:r2,3), 'o','MarkerFaceColor','b','MarkerEdgeColor','b');
            hold on
            x = linspace(d0(r1,1)-t0, d0(r2,1)-t0);
            y = zeros(size(x));
            y(:) = bkg(i,2);
            line(x./60,y,'Color','red','LineStyle','-', 'LineWidth', 2);
            xlabel("time, min");
            ylabel("gas, ppm %");
        end
        subplot(n_days, day_periods, i+1);
        plot(d0(time_pos(end):end,1)./60-t0./60, d0(time_pos(end):end,3), 'o','MarkerFaceColor','b','MarkerEdgeColor','b');
        hold on
        x = linspace(d0(time_pos(end),1)-t0, d0(end,1)-t0);
        y = zeros(size(x));
        y(:) = bkg(i+1,2);
        line(x./60,y,'Color','red','LineStyle','-', 'LineWidth', 2);
        xlabel("time, min");
        ylabel("gas, ppm %");
    end

    function plot_bkg_uncorr_obs()
        figure(4);
        t0 = min( d0(1,1), d1(1,1) );
        r2 = 0;
        for i = 1:size(pos_arr,1)-1
            r1 = r2 + 1;
            r2 = time_pos(i);
            subplot(n_days, day_periods, i);
            hold on
            plot(d0(r1:r2,1)./60-t0./60, d0(r1:r2,3), 'o','MarkerFaceColor','b','MarkerEdgeColor','b');
            plot(d1(pos_arr(i):pos_arr(i+1),1)./60-t0./60, d1(pos_arr(i):pos_arr(i+1),3), 'o','MarkerFaceColor','r','MarkerEdgeColor','r');
            
            x = linspace(d0(r1,1)-t0, d0(r2,1)-t0);
            y = zeros(size(x));
            y(:) = bkg(i,2);        
            line(x./60,y,'Color','k','LineStyle','-', 'LineWidth', 2);
            
            xlabel("time, min");
            ylabel("gas, ppm %");
        end
    end

    function bkg_arr = get_background( arr, pos_arr, num_wnds )
        t_col = 1;
        g_col = 3;
        bkg_arr = zeros(num_wnds,2);
        r22 = 0;
        for ii = 1:num_wnds-1
            r11 = r22 + 1;
            r22 = pos_arr(ii);
            bkg_arr(ii,2) = get_mode( arr(r11:r22,g_col) ); % mode
            bkg_arr(ii,1) = arr(r11,t_col); % time
        end
        bkg_arr(end,2) = get_mode( arr(pos_arr(end):end,g_col) );
        bkg_arr(end,1) = arr(pos_arr(end),t_col);
    end

    function trait = get_trait( arr, skip, method )
        n_rec = size(arr,1);
        to_skip = ceil(n_rec * skip);
        y = arr(to_skip:end-to_skip,3); % observations
        x = arr(to_skip:end-to_skip,1); % time

        % if arr(1,2) == 2471
        %     arr(1,2)
        % end   
        if numel(y) < 10 % minimal number of records
            trait = 0;
            return;
        end
        if isempty(y) % probably there are not enough observvations
            disp(["There are not enough observations for id ", num2str(arr(1,2))]);
            y = arr(:,3);
            if numel(y) < 10 % minimal number of records
                trait = 0;
                return;
            end
        end

        duration_min = ( x(end) - x(1) )./60; % duration of record in min

        switch method
            case "integral"
                trait = trapz(x,y)/duration_min;
            case "peaks"
                % need to avoid duplicated time
                xdif = x(2:end)-x(1:end-1);
                xmask = xdif > 0;
                x = x(xmask);
                y = y(xmask);

                pks = findpeaks(y, x, "MinPeakProminence", 0.0);
                trait = sum(pks)/duration_min;
            case "peaksintegral"
                % need to avoid duplicated time
                xdif = x(2:end)-x(1:end-1);
                xmask = xdif > 0;
                x = x(xmask);
                y = y(xmask);
                
                [pks,~,w,~] = findpeaks(y, x, "MinPeakProminence", 0.0);
                trait = sum(pks .* w .* 0.5)/duration_min;
            case "average"
                trait = mean(y)/duration_min;
            otherwise
                disp("Incorect method for trait calculation!");
                return;
        end
    end

    function mode_value = get_mode( arr )
        [n_counts, bins] = histcounts(arr, 'BinMethod','fd'); % other BinMethods options (methods): fd, scott, sturges, sqrt, auto
        [~, index] = max(n_counts);
        mode_value = (bins(index) + bins(index + 1))/2;
    end
    
    function pos_arr = find_pos_by_value2( what_arr, where_arr )
        % The aim is to find first occuarence (indexes)
        % of 'what_arr' in the array 'where_arr'
    
        pos_arr = zeros(size(what_arr));
        s = size(what_arr,1);
        for j = 1:s
            p = find ( where_arr >= what_arr(j), 1);
            if isempty(p)
                break;
            end
            pos_arr(j) = p;
        end
    end

    function [ pos ] = find_pos_by_value( arr, len )
        % The aim is to find the consecutive positions (indexes)
        % in the array 'arr' coresponded to equal split of
        % the array by the interval of the length 'len';
        % we start from the first position of arr (pos = 1),
        % than looking for the index, arr(pos), which correspond to the value
        % next = arr(1) + len; continue until the end of the array reached.
    
        l = 1;
        p = 1;
        next = arr(1);
        t_pos = zeros(size(arr));
        s = size(arr,1);
        while p <= s
            next = next + len;
            p = find ( arr >= next, 1);
            if isempty(p)
                break;
            end
            t_pos(l) = p;
            l = l + 1;
        end
        pos = t_pos(t_pos > 0);
    end

    function id = read_signals_ids(filename, startRow, endRow)
    
        % Initialize variables.
        if nargin<=2
            startRow = 1;
            endRow = inf;
        end
        
        % Format for each line of text:
        %   column1: double (%f)
        % For more information, see the TEXTSCAN documentation.
        formatSpec = '%10f%[^\n\r]';
        
        fileID = fopen(filename,'r');

        if fileID == -1
            disp( strcat("Cannot open the file: ", filename) );
            return;
        end
        
        % Read columns of data according to the format.
        % This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
        dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for block=2:length(startRow)
            frewind(fileID);
            dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            dataArray{1} = [dataArray{1};dataArrayBlock{1}];
        end
        
        fclose(fileID);
        
        id = [dataArray{1:end-1}];
    
    end

    function out = filter_signals( in_arr, list )
        out = in_arr;
        smask = out(:, end) == list(1,1);
        for k = 2:size(list,1)
            smask2 = out(:, end) == list(k,1);
            smask = smask + smask2;
        end
        smask = logical(table2array(smask));        
        out(~smask, :) = [];
    end


end % end of the parent function