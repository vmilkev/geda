function make_report( this, type, text, num_dat, map_dat )

if strcmp(type,"txt")
    report1;
elseif strcmp(type,"skew")
    report7;
elseif strcmp(type,"stats")
    report8;
elseif strcmp(type,"dat")
    report3;
elseif strcmp(type,"datint")
    report3int;
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

    function report3
        fprintf(this.f_hdl,'%s\n','');
        if ( numel(num_dat) == 2 )
            fprintf(this.f_hdl,'%s %12d %12d\n',text, uint64(num_dat));
        else
            %fprintf(this.f_hdl,'%s %d\n',text, uint64(num_dat));
            fprintf(this.f_hdl,'%s %.3f\n',text, (num_dat));
        end
    end

    function report3int
        fprintf(this.f_hdl,'%s\n','');
        if ( numel(num_dat) == 2 )
            fprintf(this.f_hdl,'%s %12d %12d\n',text, uint64(num_dat));
        else
            fprintf(this.f_hdl,'%s %d\n',text, uint64(num_dat));
        end
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