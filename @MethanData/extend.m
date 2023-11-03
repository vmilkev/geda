function [ d3, d4 ] = extend( this, d, d_ns )

    % Input data: d( serial_time, co2 ).
    
    tosec = 24*60*60; % constant
    
    % calculate number of samples needed to cover the
    % range_sec = ( d(end,1)-d(1,1) )*tosec, assuming (smallest)
    % 1 sec sampling interval
    samples = int64( ( d(end,1)-d(1,1) )*tosec );

    if (samples <= 0)
        disp('extend(): the samples variable is negative or equals 0.');
    end

    this.deltaT_original = median((d(2:end,1)-d(1:end-1,1))*tosec);
    
    if (this.deltaT_original <= 0)
        disp('extend(): the this.deltaT_original variable is negative or equals 0.');
    end

    d3 = zeros( samples,2 ); % extended, based on d(:,:), time series data (output)
    d4 = zeros( samples,size(d_ns,2) ); % extended, based on d(:,:), time series data (output)
    
    for i = 1:size(d3,1)
        % assign and transform time in seconds to a serial time,
        % sampling interval is 1 sec.
        d3(i,1) = ( (i-1)/tosec ) + d(1,1);
        d4(i,1) = d3(i,1);
    end
    
    % copy TRUE d values to the relevant positions in d3,
    % here we use (required) serial time representation.
    for i = 1:size(d,1)

        t = ( d(i,1)-d(1,1) )*tosec;
        j = int64(t) + 1;
        
        d3(j,1) = ( t/tosec ) + d(1,1); % time
        d4(j,1) = d3(j,1);              % time

        d3(j,2) = d(i,2);               % co2
        
%         d4(j,2) = d_ns(i,2);            % ch4
%         d4(j,3) = d_ns(i,3);            % co2
        d4(j,2:end) = d_ns(i,2:end);
    end
    
    %return;
    % ------------------------------------------------------------------
    % Time-value interpolation
    
    minval = min( d(:,2) )*0.95;

%     minval_co2 = min( d_ns(:,3) )*0.95;
%     minval_ch4 = min( d_ns(:,2) )*0.95;
    minval_d4_all = min( d_ns(:,2:end) ).*0.95;
    
    % The first loop is provide interpolation for a relatively small gaps,
    % where empty intervals between the true records less then 5 sec.
    i = 1;
    while i < size(d3,1)
        v = d3(i,2);
%         v_co2 = d4(i,3);
%         v_ch4 = d4(i,2);
        v_d4_all = d4(i,2:end);

        nzel = 0;
        if v ~= 0.0
            i1 = i;
            i = i + 1;
            while d3(i,2) == 0.0 && nzel <= 5
                nzel = nzel + 1;
                i = i + 1;
            end
            v2 = d3(i,2);
%             v2_co2 = d4(i,3);
%             v2_ch4 = d4(i,2);
            v2_d4_all = d4(i,2:end);
            i2 = i;
            if isequal(sign(v), sign(v2)) && nzel > 0
                d3(i1+1:i2-1,2) = (v+v2)/2;
%                 d4(i1+1:i2-1,2) = (v_ch4+v2_ch4)/2;
%                 d4(i1+1:i2-1,3) = (v_co2+v2_co2)/2;
                for j1 = i1+1:i2-1
                    d4(j1,2:end) = (v_d4_all+v2_d4_all)./2;
                end
            end
        else
            i = i + 1;
        end
    end
    
    % The second loop covers the rest gaps;
    % if a gap is longer then 200 sec, interpolated values are equal 'minval'.
    i = 1;
    while i < size(d3,1)
        v = d3(i,2);
%         v_co2 = d4(i,3);
%         v_ch4 = d4(i,2);
        v_d4_all = d4(i,2:end);

        nzel = 0;
        if v ~= 0.0
            i1 = i;
            i = i + 1;
            while d3(i,2) == 0.0
                nzel = nzel + 1;
                i = i + 1;
            end
            v2 = d3(i,2);
%             v2_co2 = d4(i,3);
%             v2_ch4 = d4(i,2);
            v2_d4_all = d4(i,2:end);
            i2 = i;
            if isequal(sign(v), sign(v2)) && nzel > 0
                if nzel > 200
                    d3(i1+1:i2-1,2) = minval;
%                     d4(i1+1:i2-1,2) = minval_ch4;
%                     d4(i1+1:i2-1,3) = minval_co2;
                    for j1 = i1+1:i2-1
                        d4(j1,2:end) = minval_d4_all;
                    end
                else
                    d3(i1+1:i2-1,2) = (v+v2)/2;
%                     d4(i1+1:i2-1,2) = (v_ch4+v2_ch4)/2;
%                     d4(i1+1:i2-1,3) = (v_co2+v2_co2)/2;
                    for j1 = i1+1:i2-1
                        d4(j1,2:end) = (v_d4_all+v2_d4_all)./2;
                    end
                end
            end
        else
            i = i + 1;
        end
    end

    for j2 = 2:samples
        if (d4(j2,2) == 0)
            d4(j2,2:end) = d4(j2-1,2:end);
        end
        if (d3(j2,2) == 0)
            d3(j2,2) = d3(j2-1,2);
        end
    end

end