function [A, A2, range] = resample_ams( this, data )
% 
% Transform (reshape) AMS data to a time series with resolution delta_t
%
% INPUT:
%
% data - data to transform to TS
%       data(:,1) - ID
%       data(:,2) - T_in
%       data(:,3) - T_out
%
% delta_t - time interval between records, sec
%
% OUTPUT:
%
% A(:,1) - time
% A(:,2) - value {0,1}
% A(:,3) - ID
%
% A2 is the same as A but resampled at 1sec frequency

delta_t = 1; % resample data for 1 sec first

if delta_t == 0
    this.make_report("dat", "WARNING: ERROR in create_ts2(...): delta_t == 0!", []);
    return;
end

delta_t = delta_t / (24 * 60 * 60); % transform seconds to serial seconds

c_id = 1;
t_in = 2; % assume T_in comes first
t_out = 3; % assume T_out comes last

sz = size(data);

range = data(end,t_out) - data(1,t_in);

if (range <= 0)
    this.make_report("dat", "WARNING: AMS data range (between start and end of records) is negative or equals 0!", []);
    return;
end

max_el = floor(range*2.0/delta_t);
A = zeros( max_el,3 );

loops = 0;
for i = 1:sz(1,1)
    % data for specific Cow_ID, will be associated with 1.0
    Trange = data(i,t_out) - data(i,t_in); % get range of the record

    jloops = floor(Trange/delta_t); % number of loops required to to have 'delta' sec step between the records
    loops = loops + 1;
    A(loops,3) = data(i,c_id); % last column is a cow_ID record
    A(loops,1) = data(i,t_in); % T_in, the first T record for cow_ID
    A(loops,2) = 1.0;
    for j = 1:jloops
        loops = loops + 1;
        A(loops,3) = data(i,c_id);
        A(loops,1) = A(loops-1,1) + delta_t;
        A(loops,2) = 1.0;
    end
    
    if ( i == sz(1,1) )
        break;
    end
    
    % data between specific Cow_IDs records, will be associated with 0.0
    Trange = data(i+1,t_in) - data(i,t_out); % get range of the record

    jloops = floor(Trange/delta_t);
    for j = 1:jloops
        loops = loops + 1;
        A(loops,3) = 0;
        A(loops,1) = A(loops-1,1) + delta_t;
        A(loops,2) = 0.0;
    end
end

if (loops < max_el)
    A = A(1:loops,:);
end

A2 = A; % 1 sec frequency data

delta_t = this.sampl_res_const;

if delta_t > 1 % because we already have the data resampled for 1 sec frequency
    % create resampling indeces
    l = 1;
    for i = 1:delta_t:size(A2,1)
        samples(l,1) = uint32(i);
        l = l + 1;
    end
    % resampling
    A = A2(samples,:); % only for A because we will need A2 only with 1 sec frequency
end
