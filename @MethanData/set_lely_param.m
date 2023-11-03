function set_lely_param( this, in_param )

% Get necessary properties of LelyData file and set lely_param structure
% This is overloaded method. The IN argument can be
% 1) string/char array (file name) or 2) structure of lely_param type.

if ( isstring(in_param) || ischar(in_param) ) % check if IN argument is a file name
    this.lely_param = this.read_fparam( in_param );
elseif ( isstruct(in_param) ) % otherwise, check if it is a structure
    this.lely_param = in_param;
else
    disp("There is some problem with supplied parameter's TYPE/VALUE.");
    return;
end

% Are there enough parameters in lely_param?
if ( isempty(this.lely_param.amsfile) )
    disp("The file with AMS data is not specified!");
    return;
elseif ( isempty(this.lely_param.amsrange) )
    disp("The range for AMS data is not specified!");
    return;
elseif ( isempty(this.lely_param.device) )
    disp("The device name (robot) of AMS data is not specified!");
    return;
elseif ( isempty(this.lely_param.id) )
    disp("The animal ID name of AMS data is not specified!");
    return;
elseif ( isempty(this.lely_param.time) )
    disp("The TIME option for AMS data is not specified!");
    return;
end
    
end
