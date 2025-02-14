function [ d3, d4, data_range_hours, sampl_frq_median, sampl_frq_min, sampl_frq_max ] = resample_snf( this, d, d_ns )

    % Input data: d( serial_time, co2 ).
    
    tosec = 24*60*60; % constant
    
    % calculate number of samples needed to cover the
    % range = ( d(end,1)-d(1,1) )*tosec, assuming a 1 sec sampling interval
    
    samples = int64( ( d(end,1)-d(1,1) )*tosec ); % the range

    if (samples <= 0)
        this.make_report("dat", "WARNING in extend(): The sniffers data range is negative or equals 0!", []);
        return;
    end

    this.deltaT_original = median((d(2:end,1)-d(1:end-1,1))*tosec); % actual sampling frequency of the data
    
    if (this.deltaT_original <= 0)
        this.make_report("dat", "WARNING in extend(): The median sampling frequency of sniffers data is negative or equals 0!", []);
        return;
    end

    data_range_hours = samples/3600;
    sampl_frq_median = this.deltaT_original;
    sampl_frq_min = min((d(2:end,1)-d(1:end-1,1))*tosec);
    sampl_frq_max = max((d(2:end,1)-d(1:end-1,1))*tosec);
    
    % extend the entire data range using 1 sec sampling interval;
    % arrays d and d_ns will be mapped into d3 and d4 with d3(i,:) == 0
    % if there is no specific recors in d of d_ns which can be mapped
    % into d3(i,:); this will be resampled (interpolated)
    
    d3 = zeros( samples,2 ); % extended, based on d(:,:), time series data (output)
    d4 = zeros( samples,size(d_ns,2) ); % extended, based on d(:,:), time series data (output)
        
    % copy TRUE d values to the relevant positions in d3,
    % here we use (required) serial time representation.

    for i = 1:size(d,1)
        t = ( d(i,1)-d(1,1) )*tosec;
        j = int64(t) + 1;
        d3(j,1) = ( t/tosec ) + d(1,1); % convert and assign time in seconds to a serial time,
        d4(j,1) = d3(j,1);              % time
        d3(j,2) = d(i,2);               % co2        
        d4(j,2:end) = d_ns(i,2:end);    % co2 & ch4
    end

    non_zerros = find( d3(:,1) ~= 0 ); % find missing records

    if non_zerros(1) ~= 1 % if the very first records are missing
        for i = 1:non_zerros(1)-1
            d3(i,1) = ( (i-1)/tosec ) + d(1,1); % time
            d4(i,1) = d3(i,1); % time
            d3(i,2) = d3(non_zerros(1),2); % co2 
            d4(i,2:end) = d4(non_zerros(1),2:end); % co2 & ch4
        end
    end

    for i = 1:size(non_zerros,1)-1 % interpolate data at missing locations
        if ( non_zerros(i+1) - non_zerros(i) ) > 1
            x_1 = d3(non_zerros(i),1);
            x_2 = d3(non_zerros(i+1),1);
            y3_1 = d3(non_zerros(i),2);
            y3_2 = d3(non_zerros(i+1),2);
            y4_1 = d4(non_zerros(i),2:end);
            y4_2 = d4(non_zerros(i+1),2:end);
            for j = non_zerros(i)+1:non_zerros(i+1)-1
                d3(j,1) = ( (j-1)/tosec ) + d(1,1); % time
                d4(j,1) = d3(j,1); % time
                d3(j,2) = interp_missing( y3_1, y3_2, x_1, x_2, d3(j,1) ); % co2 
                d4(j,2:end) = interp_missing( y4_1, y4_2, x_1, x_2, d3(j,1) ); % co2 & ch4
            end
        end
    end
    
    clear non_zerros

    % here we define which sampling rate to be applied
    use_sampl_rate = this.sampl_res_const;
    if use_sampl_rate == 0 % if the rate is not set (default is 0)
        use_sampl_rate = this.deltaT_original; % use the data actual rate
        this.sampl_res_const = this.deltaT_original; % this is needed for AMS data
    end

    if use_sampl_rate > 1 % because we already have the data resampled for 1 sec frequency
        % create resampling indeces
        l = 1;
        for i = 1:use_sampl_rate:size(d3,1)
            samples(l,1) = uint32(i);
            l = l + 1;
        end
        % resampling
        d3 = d3(samples,:); % only for d3 because we will need d4 only with 1 sec frequency
        %d4 = d4(samples,:);
    end
    
    % ------ INTERNAL FUNCTION ------------------

    function y = interp_missing(y1,y2,x1,x2,x)
        y = (x-x1)*(y2-y1)/(x2-x1) + y1;
    end

end