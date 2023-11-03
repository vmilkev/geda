function set_mesr_param( this, in_param )

% Get necessary properties of MesrData file and set mesr_param structure
% This is overloaded method. The IN argument can be
% 1) string/char array (file name) or 2) structure of mesr_param type.

if ( isstring(in_param) || ischar(in_param) ) % check if IN argument is a file name
    this.mesr_param = this.read_fparam( in_param );
elseif ( isstruct(in_param) ) % otherwise, check if it is a structure
    this.mesr_param = in_param;
else
    disp("There is some problem with supplied parameter's TYPE/VALUE.");
    return;
end

% Are there enough parameters in lely_param?
if ( isempty(this.mesr_param.sniffile) )
    disp("The file with Sniffer data is not specified!");
    return;
elseif ( isempty(this.mesr_param.snifrange) )
    disp("The range for Sniffer data is not specified!");
    return;
elseif ( isempty(this.mesr_param.robid) )
    disp("The device name (robot) of Sniffer data is not specified!");
    return;
end

this.device = this.mesr_param.robid;

end
