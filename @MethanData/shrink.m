function B = shrink( this, data, delta_t )

% function to reduce time interval between records
% 
% INPUT:
% 
% data - data file
%        data(:,1) - time
%        data(:,2) - record CO2
%        data(:,3) - record CH4 <- this one is not present in scaled data!
%        
% delta - required time interval, sec
% 
% OUTPUT:
% 
% B - shrinked data

if ~exist('delta_t','var') || isempty(delta_t)
    delta_t = 0.00005787; % 5 sec
else
    delta_t = delta_t / (24 * 60 * 60); % transform seconds to serial seconds
end

% get records interval
szB = size(data);

dt = mean( data(2:end,1) - data(1:end-1,1) );

interval = ( uint64( delta_t/dt ) );

index =  ~mod( uint64( (1:szB(1,1)) )',interval ) ;

B = data(index,:);

end