function make_phenotypes( fname, w_length )

close all

% Takes as an input the GEDA produced file
% (such as aligned_init_sniffer_101.txt)
% and do the following (for each gas type present):
% 1. calculates the backgroung values;
% 2. correct data for the background
% 3. do SSA-based denoising
% 4. calculates emission phenotype.
%
% fname - input file name
% w_length - max duration in hours used for background calculations

if nargin < 1
    disp("At least a GEDA resulting file with aligned data is expected as an input argument!");
    return;
end

if nargin < 2
    
    if ~exist('w_length','var')
      w_length = 6;
    end
    
    if ~exist('fname','var')
        disp("At least a GEDA resulting file with aligned data is expected as an input argument!");
        return;
    end
end

[res] = read_aligned(fname);

time_col = 1;
id_col = 2;
outlier_low_limit = -0.1;

% determine the number of gas types
res_dim = size(res);
n_gases = res_dim(1,2) - 3; % 3 because: time, id, and signal cols
max_data_rows = res_dim(1,1);

for igas = 3:3%res_dim(1,2)-1 % loop over gas types

    range = 1:floor(max_data_rows/10); % for testing/debugging use a part of the data
    %range = 1:max_data_rows; % use the entire range
    
    % convert tables to arrays
    gas = table2array(res(range,igas));
    id = table2array(res(range,id_col));
    t = convertTo(table2array(res(range,1)), 'posixtime');
    % ------------------------
    
    mask = find(gas >= outlier_low_limit); % mask outliers
        
    d = [t(mask) id(mask) gas(mask)]; % all data without filtered outliers
    
    mask = d(:,id_col) == 0; % mask the background dataset
    
    d0 = d(mask,:); % background dataset
    d1 = d(~mask,:); % not-background dataset

    clear d
    
    day_in_sec = 24*60*60; % duration of one day, seconds
    n_days = ceil( ( d0(end, 1) - d0(1,1) ) / day_in_sec ); % number of days in the dataset, days
    wnd_len = min(w_length*60*60, d0(end, 1) - d0(1,1)); % duration of moving window, seconds
    
    time_pos  = find_pos_by_value( d0(:,1), wnd_len );
    
    n_wnd = size(time_pos,1)+1; % num of windows for the entire record
    day_periods = ceil(n_wnd/n_days); % num of windows per day
       
    bkg = get_background( d0, time_pos, n_wnd ); % calculate background concentration for each window
    pos_arr = find_pos_by_value2( bkg(:,1), d1(:,1) );
    
    d = d1; % make copy in order to have a plot of uncorrected d1
    for i = 1:size(pos_arr,1)-1
        d(pos_arr(i):pos_arr(i+1)-1,3) = d(pos_arr(i):pos_arr(i+1)-1,3) - bkg(i,2);
    end
    d(pos_arr(end):end,3) = d(pos_arr(end):end,3) - bkg(end,2);

    mask = d(:,3) > 0; % find negative outliers
    d = d(mask,:); % exclude negative outliers

    unique_ids = unique( d(:,2) );
    for i = 1:size(unique_ids,1) % loop over specific ids
        mask = d(:,2) == unique_ids(i);
        d_id = d(mask,:);
        t_gaps = d_id(2:end,1) - d_id(1:end-1,1);
        md = get_mode(t_gaps);
        gaps = find( t_gaps > md * 10 );
        if ~isempty(gaps)
            r1 = 1;
            for h = 1:size(gaps,1) % loop over distant records of a specific id
                i_range = r1:gaps(h);
                if h == size(gaps,1)
                    i_range = gaps(h) + 1:size(d_id,1);
                end
                record = d_id(i_range, :);
                r1 = gaps(h) + 1;
            end
        else % only one record is present
            record = d_id;
        end
    end


    % figure(1);
    % t0 = d0(1,1);
    % plot(bkg(:,1)./60 -t0./60, bkg(:,2), '-o');
    
    figure(2)
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

    figure(3);
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


end % end of for_loop ove gas types

% NESTED FUNCTIIONS

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

end % end of the parent function