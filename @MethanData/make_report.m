function make_report( this, type, text, num_dat, map_dat, num_dat2, map_dat2 )

if strcmp(type,"txt")
    report1;
elseif strcmp(type,"skew")
    report7;
elseif strcmp(type,"stats")
    report8;
elseif strcmp(type,"fig")
    if nargin == 4
        report2;
                 % num_dat == msr_num
    elseif nargin == 5
        report4;
                 % num_dat == msr_ts
                 % map_dat == lely_ts
    elseif nargin > 5
        report5; % plot all robot info & additional mesurements info
                 % num_dat == mesr_num,
                 % map_dat == mesr_map,
                 % num_dat2 == lely_num,
                 % map_dat2 == lely_map 
    end
elseif strcmp(type,"dat")
    report3;
elseif strcmp(type,"res")
    report6;
end

%% FUNCTIONS:

    function report1
        disp_rows = 5; % number of rows to print
        fprintf(this.f_hdl,'%s\n',text);
        fprintf(this.f_hdl,'%s\n','');
        fprintf(this.f_hdl,'%s %d %s\n','First and last', disp_rows, 'rows in data file:');
        fprintf(this.f_hdl,'%s\n','');
        if ( numel(map_dat) > 5 )
            fprintf(this.f_hdl,'%s\n','(... only the first 5 datasets will be displayed.)');
        end
        for i = 1:numel(map_dat)
            if (i > 5) % (... only the first 5 datasets will be displayed.)
                break;
            end
            max_el = size(map_dat{i},1);
            if (max_el > 2*disp_rows)
                fprintf(this.f_hdl,'%s %d\n','Table data:',i);
                writetable(map_dat{i}(1:disp_rows,:), this.text_fname, 'Delimiter', 'tab', 'WriteMode', 'append');
                fprintf(this.f_hdl,'%s\n','');
                writetable(map_dat{i}(end-disp_rows:end,:), this.text_fname, 'Delimiter', 'tab', 'WriteMode', 'append');
                fprintf(this.f_hdl,'%s\n','');
            else
                fprintf(this.f_hdl,'%s %d\n','Table data:',i);
                writetable(map_dat{i}(1:end,:), this.text_fname, 'Delimiter', 'tab', 'WriteMode', 'append');
                fprintf(this.f_hdl,'%s\n','');
            end
        end
        fprintf(this.f_hdl,'%s\n','----------------------------------------------------------------');
        fprintf(this.f_hdl,'%s\n','');
        fprintf(this.f_hdl,'%s %d %s\n','First and last', disp_rows, 'rows of processed numeric data from file:');
        fprintf(this.f_hdl,'%s\n','');
        for i = 1:numel(map_dat)
            if (i > 5) % (... only the first 5 datasets will be displayed.)
                break;
            end
            max_el = size(map_dat{i},1);
            if (max_el > 2*disp_rows)
                fprintf(this.f_hdl,'%s %d\n','Numeric data:',i);
                writematrix(num_dat{i}(1:disp_rows,:), this.text_fname, 'Delimiter', 'tab', 'WriteMode', 'append');
                fprintf(this.f_hdl,'%s\n','');
                writematrix(num_dat{i}(end-disp_rows:end,:), this.text_fname, 'Delimiter', 'tab', 'WriteMode', 'append');
                fprintf(this.f_hdl,'%s\n','');
            else
                fprintf(this.f_hdl,'%s %d\n','Numeric data:',i);
                writematrix(num_dat{i}(1:end,:), this.text_fname, 'Delimiter', 'tab', 'WriteMode', 'append');
                fprintf(this.f_hdl,'%s\n','');
            end
        end
        fprintf(this.f_hdl,'%s\n','----------------------------------------------------------------');
        fprintf(this.f_hdl,'%s\n','');        
    end

    function report2
        f1 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(num_dat)
            plot( num_dat{i,1}(:,1), num_dat{i,1}(:,2) );
        end
        hold off;
        box on;
        title(['CH4 vs Time; robot no. ', num2str(this.get_mesr_robot())]);
        xlabel('serial time');
        ylabel('CH4');

        f2 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(num_dat)
            histogram(num_dat{i,1}(:,2),'Normalization','probability');
        end
        hold off;
        box on;
        title(['CH4 distribution; robot no. ', num2str(this.get_mesr_robot())]);
        ylabel('frequency');
        xlabel('CH4');

        f3 = figure('visible','off');
        set(gcf,'color','w');
        %
        subplot(2,2,1);
        hold on;
        for i = 1:numel(num_dat)    
            plot( num_dat{i,1}(1:end-1,2), num_dat{i,1}(2:end,2), '.' );
        end
        hold off;
        box on;
        title('CH4 lags plot');
        xlabel('$CH{4}(T{i})$','interpreter','latex');
        ylabel('$CH{4}(T{i+1})$','interpreter','latex');
        %
        subplot(2,2,2);
        hold on;
        for i = 1:numel(num_dat)    
            plot( num_dat{i,1}(1:end-2,2), num_dat{i,1}(3:end,2), '.' );
        end
        hold off;
        box on;
        title('CH4 lags plot');
        xlabel('$CH{4}(T{i})$','interpreter','latex');
        ylabel('$CH{4}(T{i+2})$','interpreter','latex');
        %
        subplot(2,2,3);
        hold on;
        for i = 1:numel(num_dat)    
            plot( num_dat{i,1}(1:end-3,2), num_dat{i,1}(4:end,2), '.' );
        end
        hold off;
        box on;
        title('CH4 lags plot');
        xlabel('$CH{4}(T{i})$','interpreter','latex');
        ylabel('$CH{4}(T{i+3})$','interpreter','latex');
        %
        subplot(2,2,4);
        hold on;
        for i = 1:numel(num_dat)    
            plot( num_dat{i,1}(1:end-4,2), num_dat{i,1}(5:end,2), '.' );
        end
        hold off;
        box on;
        title('CH4 lags plot');
        xlabel('$CH{4}(T{i})$','interpreter','latex');
        ylabel('$CH{4}(T{i+4})$','interpreter','latex');

        % ---------------------------------------------------
        f4 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(num_dat)
            plot( num_dat{i,1}(:,1), num_dat{i,1}(:,3) );
        end
        hold off;
        box on;
        title(['CO2 vs Time; robot no. ', num2str(this.get_mesr_robot())]);
        xlabel('serial time');
        ylabel('CO2');

        f5 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(num_dat)
            histogram(num_dat{i,1}(:,3),'Normalization','probability');
        end
        hold off;
        box on;
        title(['CO2 distribution; robot no. ', num2str(this.get_mesr_robot())]);
        ylabel('frequency');
        xlabel('CO2');

        f6 = figure('visible','off');
        set(gcf,'color','w');
        %
        subplot(2,2,1);
        hold on;
        for i = 1:numel(num_dat)    
            plot( num_dat{i,1}(1:end-1,3), num_dat{i,1}(2:end,3), '.' );
        end
        hold off;
        box on;
        title('CO2 lags plot');
        xlabel('$CO{2}(T{i})$','interpreter','latex');
        ylabel('$CO{2}(T{i+1})$','interpreter','latex');
        %
        subplot(2,2,2);
        hold on;
        for i = 1:numel(num_dat)    
            plot( num_dat{i,1}(1:end-2,3), num_dat{i,1}(3:end,3), '.' );
        end
        hold off;
        box on;
        title('CO2 lags plot');
        xlabel('$CO{2}(T{i})$','interpreter','latex');
        ylabel('$CO{2}(T{i+2})$','interpreter','latex');
        %
        subplot(2,2,3);
        hold on;
        for i = 1:numel(num_dat)    
            plot( num_dat{i,1}(1:end-3,3), num_dat{i,1}(4:end,3), '.' );
        end
        hold off;
        box on;
        title('CO2 lags plot');
        xlabel('$CO{2}(T{i})$','interpreter','latex');
        ylabel('$CO{2}(T{i+3})$','interpreter','latex');
        %
        subplot(2,2,4);
        hold on;
        for i = 1:numel(num_dat)    
            plot( num_dat{i,1}(1:end-4,3), num_dat{i,1}(5:end,3), '.' );
        end
        hold off;
        box on;
        title('CO2 lags plot');
        xlabel('$CO{2}(T{i})$','interpreter','latex');
        ylabel('$CO{2}(T{i+4})$','interpreter','latex');

        print(f1, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f2, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f3, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f4, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f5, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f6, '-dpsc', this.fig_resl, '-append', this.figs_fname);
    end

    function report3
        fprintf(this.f_hdl,'%s\n','');
        if ( numel(num_dat) == 2 )
            fprintf(this.f_hdl,'%s %12d %12d\n',text, uint64(num_dat));
        else
            %fprintf(this.f_hdl,'%s %d\n',text, uint64(num_dat));
            fprintf(this.f_hdl,'%s %.3f\n',text, (num_dat));
        end
    end

    function report4
        % num_dat == msr_ts
        % map_dat == lely_ts
        
        f1 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(num_dat)
            plot(num_dat{i, 1}(:,1),num_dat{i, 1}(:,2), '-');
        end
        plot(map_dat{1, 1}(:,1), map_dat{1, 1}(:,2), '.');
        hold off;
        box off;
        %ylim([-0.1 1.1]);
        l = string([1:numel(num_dat)]');
        l2 = "ams data";
        legend([l;l2]);
        title(['Binary AMS & normalized CO2 signals; robot no. ', num2str(this.get_mesr_robot())]);
        xlabel('$T$','interpreter','latex');
        ylabel('$Signal$','interpreter','latex');

        f2 = figure('visible','off');
        set(gcf,'color','w');
        subplot(2,1,1);
        histogram(map_dat{1,1}(:,2),'Normalization','probability');
        title(['Frequencys of binary AMS signal ; robot no. ', num2str(this.get_mesr_robot())]);
        ylabel('Frequency');
        xlabel('values of AMS signal');

        subplot(2,1,2);
        hold on;
        for i = 1:numel(num_dat)
            histogram(num_dat{i,1}(:,2),'Normalization','probability');
        end
        hold off;
        box on;
        title(['Frequencys of normalized CO2 signal; robot no. ', num2str(this.get_mesr_robot())]);
        ylabel('Frequency');
        xlabel('values of normalized CO2 signal');
        
        print(f1, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f2, '-dpsc', this.fig_resl, '-append', this.figs_fname);
    end

    function report5
         % num_dat == mesr_num,
         % map_dat == mesr_map,
         % num_dat2 == lely_num,
         % map_dat2 == lely_map 

        f2 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(this.Robots)
            if ismember('Begintime',map_dat2{i, 1}.Properties.VariableNames)
                if ~isa(map_dat2{i, 1}.Begintime,'cell')
                    histogram(map_dat2{i, 1}.Begintime);
                end
            elseif ismember('DviStartTime',map_dat2{i, 1}.Properties.VariableNames)
                if ~isa(map_dat2{i, 1}.DviStartTime,'cell')
                    histogram(map_dat2{i, 1}.DviStartTime);
                end
            end
        end
        hold off;
        box on;
        %xticks(this.Robots);
        title({'Number of AMS records for all robots'});
        ylabel('AMS records');
        xlabel('Robot name');

        f3 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(this.Robots)
            if ismember('Begintime',map_dat2{i, 1}.Properties.VariableNames)
                if ~isa(map_dat2{i, 1}.Begintime,'cell')
                    histogram(map_dat2{i, 1}.Begintime,'Normalization','probability','FaceAlpha',0.5);
                end
            elseif ismember('DviStartTime',map_dat2{i, 1}.Properties.VariableNames)
                if ~isa(map_dat2{i, 1}.DviStartTime,'cell')
                    histogram(map_dat2{i, 1}.DviStartTime,'Normalization','probability','FaceAlpha',0.5);
                end
            end
        end
        hold off;
        box on;
        legend(num2str(this.Robots));
        title({'Frequency distribution of AMS records over time for all robots'});
        ylabel('Frequency');
        xlabel('Time');

        f4 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(this.Robots)
            plot(num_dat2{i, 1}(:,3), '.');
        end
        hold off;
        box on;
        legend(num2str(this.Robots));
        title({'Number of AMS records over time'});
        ylabel('Serial time');
        xlabel('AMS records');
        
        % -----------------------------------------------------------------
        f5 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(map_dat)
            if ~isa(map_dat{i, 1}.DATE,'cell')
                plot(map_dat{i, 1}.DATE,'o');
            end
        end
        hold off;
        box on;
        title(['Sniffer records over time, robot no. ', num2str(this.get_mesr_robot())]);
        ylabel('Day');
        xlabel('Sniffer records');

        f6 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(map_dat)
            histogram(num_dat{i, 1}(:,1),'Normalization','probability','FaceAlpha',0.5);
        end
        hold off;
        box on;
        title(['Frequency distribution of sniffer records, robot no. ', num2str(this.get_mesr_robot())]);
        ylabel('Frequency');
        xlabel('Serial time');

        f7 = figure('visible','off');
        set(gcf,'color','w');
        hold on;
        for i = 1:numel(map_dat)
            plot(num_dat{i, 1}(:,1), 'o');
        end
        hold off;
        box on;
        title(['Sniffer records over time, robot no. ', num2str(this.get_mesr_robot())]);
        ylabel('Serial time');
        xlabel('Sniffer records');

        f8 = figure('visible','off');
        set(gcf,'color','w');
        histogram(num_dat2{1, 1}(:,3),'Normalization','probability','FaceAlpha',0.5);
        hold on;
        for i = 1:numel(map_dat)
            histogram(num_dat{i, 1}(:,1),'Normalization','probability','FaceAlpha',0.5);
        end
        hold off;
        box on;
        title({'Frequency distribution of AMS & sniffer records; robot no.', num2str(this.get_mesr_robot())});
        ylabel('Frequency');
        xlabel('Serial time');

%         print(f1, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f2, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f3, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f4, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f5, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f6, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f7, '-dpsc', this.fig_resl, '-append', this.figs_fname);
        print(f8, '-dpsc', this.fig_resl, '-append', this.figs_fname);
    end

    function report6
        %num_dat = skew
        %map_dat = skew_dates
         
        f = figure('visible','off');
        set(gcf,'color','w');
        %---------------------------
        datasets = size(num_dat,1);

        if size(num_dat,2) < 2
            return;
        end
        
        for which_data = 1:datasets
        
            [skewdata,outliers] = rmoutliers(num_dat(which_data,:)./60); % now in minutes
            time = map_dat(which_data,outliers == 0)./60; % now in minutes
        
            [p,S] = polyfit(time,skewdata,1);
        
            [y_fit,delta] = polyval(p,time,S);
        
            subplot(datasets,1,which_data);
        
            plot(time,skewdata,'bo', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
            hold on
            plot(time,y_fit,'r-', 'LineWidth',2);
            plot(time,y_fit+2*delta,'b--',time,y_fit-2*delta,'b--');
            s1 = 'Timing Skew ';
            s2 = strcat(s1,', sniffer data file no. ');
            s1 = strcat(s2, num2str(which_data));
            title(s1);
            % Plot residuals as lines from actual data to fitted line.
            for k = 1 : length(time)
                yActual = skewdata(k);
                yFit = y_fit(k);
                x = time(k);
                line([x, x], [yFit, yActual], 'Color', 'b', 'LineStyle', ':');
            end
        
            xlabel('Time (of sniffer records) passed from the first record in the file, min','Interpreter','latex'); % Time from the first record, min
            ylabel('Timing skew, min','Interpreter','latex'); % Timing Skew, min
            xlim([min(time) max(time)]);
            set(gca,'FontSize',12);
        end
        %---------------------------
        print(f, '-dpsc', this.fig_resl, '-append', this.figs_fname);
    end

    function report7
        sz = size(num_dat);
        fprintf(this.f_hdl,'%s\n','Starting dates of sniffer and corresponding AMS records,');
        fprintf(this.f_hdl,'%s\n','for each processed data set:');
        %fprintf(this.f_hdl,'%s\n','');
        fprintf(this.f_hdl,'%s\n','----------------------------------------------------------');
        fprintf(this.f_hdl,'%s %s %s\n','no.  |','  Sniffer Dates       |','Corresponding AMS dates');
        fprintf(this.f_hdl,'%s\n','----------------------------------------------------------');
        isignal = 1;
        for j = 1:sz(1,1)
            for i = 1:sz(1,2)
                if ( ~isempty(num_dat{j,i}) && ~isempty(map_dat{j,i}) )
                    dates_table = table(isignal, num_dat(j,i), map_dat(j,i));
                    writetable(dates_table,this.text_fname,'WriteMode','Append',...
                        'WriteVariableNames',false,'Delimiter','tab','QuoteStrings',false);
                    isignal = isignal + 1;
                end
            end
        end
        fprintf(this.f_hdl,'%s\n','----------------------------------------------------------');
        fprintf(this.f_hdl,'%s\n','');
    end

    function report8
        if ( iscell(num_dat) )
            fprintf(this.f_hdl,'%s\n','');
            fprintf(this.f_hdl,'%s\n','Results of alignment and reliability detection:');
            %fprintf(this.f_hdl,'%s\n','');            
            fprintf(this.f_hdl,'%s\n','----------------------------------------------------------------------------');
            fprintf(this.f_hdl,'%s %10s %10s %12s %10s\n','Signal no. |','Skew, sec |','Test statistic |','Test statistic 2 |', 'Reliable');
            fprintf(this.f_hdl,'%s\n','----------------------------------------------------------------------------');

            LogicalStr = {'false', 'true'};

            isignal = 1;
            for j = 1:size(num_dat{1,1},1) % loop over number of data sets in each num_dat{j}
                for i = 1:size(num_dat{1,1}{j,1},2)
                    %if ( num_dat{2}(j,i) ~= 0.0 ) %( (num_dat{2}(j,i) ~= 0.0) && (num_dat{3}(j,i) ~= 0.0) && (num_dat{4}(j,i) ~= 0.0) )
                        %fprintf(this.f_hdl,'%10d %11d %11.3f\n',isignal,num_dat{1}(j,i),num_dat{2}(j,i) );
                        %fprintf(this.f_hdl,'%6d %14d %12.3f %15.3f %16d\n',isignal,num_dat{1}(j,i),num_dat{2}(j,i),num_dat{3}(j,i),num_dat{4}(j,i) );
                        
                        fprintf(this.f_hdl,'%6d %14d %12.3f %15.3f %16s\n',isignal,num_dat{1,1}{j,1}(1,i),num_dat{2,1}{j,1}(1,i),num_dat{3,1}{j,1}(1,i),LogicalStr{num_dat{4,1}{j,1}(1,i)+1} );
                        %fprintf(this.f_hdl,'%6d %14d %12.3f %15.3f %16d\n',isignal,num_dat{1,1}{j,1}(1,i),num_dat{2,1}{j,1}(1,i),num_dat{3,1}{j,1}(1,i),num_dat{4,1}{j,1}(1,i) );
                        isignal = isignal + 1;
                    %end
                end
            end
            fprintf(this.f_hdl,'%s\n','----------------------------------------------------------------------------');
            fprintf(this.f_hdl,'%s\n','');
        else
            if size(num_dat,2) == 4
                fprintf(this.f_hdl,'%s\n','');
                fprintf(this.f_hdl,'%s\n','AMS data summary:');
                fprintf(this.f_hdl,'%s\n','-------------------------------------------------------------------------');
                fprintf(this.f_hdl,'%2s %10s %10s %13s\n','Data set no. |','Length, hours |','Resampl. freq., sec |', 'Elapsed time, sec');
                fprintf(this.f_hdl,'%s\n','-------------------------------------------------------------------------');
                for j = 1:size(num_dat,1)
                    fprintf(this.f_hdl,'%8d %15d %15d %19.3f\n',num_dat(j,1),num_dat(j,2),num_dat(j,3),num_dat(j,4) );
                end
                fprintf(this.f_hdl,'%s\n','-------------------------------------------------------------------------');
                fprintf(this.f_hdl,'%s\n','');
            elseif size(num_dat,2) == 5
                fprintf(this.f_hdl,'%s\n','');
                fprintf(this.f_hdl,'%s\n','Data processing summary:');
                fprintf(this.f_hdl,'%s\n','---------------------------------------------------------------------------------------');
                fprintf(this.f_hdl,'%2s %10s %10s %10s %10s\n','Data set no. |','Length, hours |','Num. signals |', 'Signal length, hour |', 'Elapsed time, sec');
                fprintf(this.f_hdl,'%s\n','---------------------------------------------------------------------------------------');
                for j = 1:size(num_dat,1)
                    fprintf(this.f_hdl,'%8d %15d %15d %15d %19.3f\n',num_dat(j,1),num_dat(j,2),num_dat(j,3),num_dat(j,4),num_dat(j,5) );
                end
                fprintf(this.f_hdl,'%s\n','---------------------------------------------------------------------------------------');
                fprintf(this.f_hdl,'%s\n','');
            elseif size(num_dat,2) == 6
                fprintf(this.f_hdl,'%s\n','');
                fprintf(this.f_hdl,'%s\n','Sniffer data summary:');
                fprintf(this.f_hdl,'%s\n','-------------------------------------------------------------------------------------------------------');
                fprintf(this.f_hdl,'%2s %10s %10s %10s %10s %10s\n','Data set no. |','Length, hours |','Median freq., sec |', 'Min.freq., sec |', 'Max.freq., sec |', 'Elapsed time, sec');
                fprintf(this.f_hdl,'%s\n','-------------------------------------------------------------------------------------------------------');
                for j = 1:size(num_dat,1)
                    fprintf(this.f_hdl,'%8d %15d %16d %16d %16d %17.3f\n',num_dat(j,1),num_dat(j,2),num_dat(j,3),num_dat(j,4),num_dat(j,5),num_dat(j,6) );
                end
                fprintf(this.f_hdl,'%s\n','-------------------------------------------------------------------------------------------------------');
                fprintf(this.f_hdl,'%s\n','');
            elseif size(num_dat,2) == 7
                fprintf(this.f_hdl,'%s\n','');
                fprintf(this.f_hdl,'%s\n','Signal length optimization summary:');
                fprintf(this.f_hdl,'%s\n','------------------------------------------------------------------------------------------------------');
                fprintf(this.f_hdl,'%2s %10s %10s %10s %10s %10s\n','Sig.len., hours |','Reliable, % |','Confidence, % |','Mean skew, sec |', 'Std skew, sec |', 'Elapsed time, sec');
                fprintf(this.f_hdl,'%s\n','------------------------------------------------------------------------------------------------------');
                for j = 1:size(num_dat,1)
                    if num_dat(j,1) == 0
                        continue;
                    end
                    fprintf(this.f_hdl,'%10d %15.3f %15.3f %15.3f %16.3f %14.3f\n',num_dat(j,1),num_dat(j,2),num_dat(j,3),num_dat(j,4),num_dat(j,5),num_dat(j,6) );
                end
                fprintf(this.f_hdl,'%s\n','------------------------------------------------------------------------------------------------------');
                fprintf(this.f_hdl,'%s\n','');
            end
        end
    end
    
end