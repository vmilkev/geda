function A = create_ts( this, data, header, delta_t )
% 
% use Lely data and reshape it to a time series
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

if ~exist('delta_t','var') || isempty(delta_t)
    delta_t = 0.00005787; % 5 sec
else
    delta_t = delta_t / (24 * 60 * 60); % transform seconds to serial seconds
end

% find indexes for ID, Tin, Tout in the 'data'
c_id = find(header == "id");
alltime = find(header == "time");
t_in = alltime(1); % assume T_in comes first
t_out = alltime(2); % assume T_out comes last

sz = size(data);

range = data(end,t_out) - data(1,t_in);

if (range <= 0)
    disp('create_ts(): the data range is negative or equals 0.');
end

max_el = floor(range*2.0/delta_t);
A = zeros( max_el,3 );

loops = 0;
for i = 1:sz(1,1)
    % data for specific Cow_ID, will be associated with 1.0
    Trange = data(i,t_out) - data(i,t_in); % get range of the record

    jloops = floor(Trange/delta_t); % number of loops required to to have 'delta' sec step between the records
%     if ( loops >= 1 && t_in < A(loops,1) )
%         loops = loops - 1;
%     end
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
%      loops = loops + 1;
%      A(loops,3) = data(i,c_id);
%      A(loops,1) = data(i,t_out); % T_out, the last T record for cow_ID
%      A(loops,2) = 1.0;
    
    if ( i == sz(1,1) )
        break;
    end
    
    % data between specific Cow_IDs records, will be associated with 0.0
    Trange = data(i+1,t_in) - data(i,t_out); % get range of the record
    %Trange = data(i+1,t_in) - A(loops,1); % get range of the record

    jloops = floor(Trange/delta_t);
    for j = 1:jloops
        loops = loops + 1;
        A(loops,3) = 0;
        A(loops,1) = A(loops-1,1) + delta_t;
        A(loops,2) = 0.0;
    end
    %disp(["loops = ", num2str(loops)]);
end

if (loops < max_el)
    A = A(1:loops,:);
end