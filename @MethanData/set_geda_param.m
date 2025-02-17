function [ sampT, siglen ] = set_geda_param( this, in_param )

% Get necessary properties of AMS file and set lely_param structure
% This is overloaded method. The IN argument can be
% 1) string/char array (file name) or 2) structure of lely_param type.

if ( isstring(in_param) || ischar(in_param) ) % check if IN argument is a file name
    this.geda_param = this.read_fparam( in_param );
elseif ( isstruct(in_param) ) % otherwise, check if it is a structure
    this.geda_param = in_param;
else
    this.make_report("dat", "WARNING: There is some problem with supplied parameter's TYPE/VALUE.", []);
    return;
end

sampT = [];
if ~isempty(this.geda_param.deltaT)
    sampT = this.geda_param.deltaT;
end

siglen = 0;
if ~isempty(this.geda_param.siglength)
    siglen = this.geda_param.siglength;
end

if isempty(this.geda_param.denoise) || this.geda_param.denoise == 0
    this.is_ssa = false; % change if not the default
else % if the ssa remains true
    if ~isempty(this.geda_param.wndlength)
        this.ssa_wnd = this.geda_param.wndlength;
    end
    if ~isempty(this.geda_param.eigenvalues)
        this.ssa_eigval = this.geda_param.eigenvalues;
    end
    if this.ssa_wnd < this.ssa_eigval
        this.make_report("dat", "WARNING: Parameter SSA1 cannot be lower than SSA2!", []);
        return;
    end
end

% These are for AMS
if ( isempty(this.geda_param.amsfile) )
    this.make_report("dat", "WARNING: The file with AMS data is not specified!", []);
    return;
elseif ( isempty(this.geda_param.amsrange) )
    this.make_report("dat", "WARNING: The range for AMS data is not specified!", []);
    return;
elseif ( isempty(this.geda_param.device) )
    this.make_report("dat", "WARNING: The device name (robot) of AMS data is not specified!", []);
    return;
elseif ( isempty(this.geda_param.id) )
    this.make_report("dat", "WARNING: The animal ID name of AMS data is not specified!", []);
    return;
elseif ( isempty(this.geda_param.time) )
    this.make_report("dat", "WARNING: The TIME option for AMS data is not specified!", []);
    return;
end

% These are for Sniffer
if ( isempty(this.geda_param.sniffile) )
    this.make_report("dat", "WARNING: The file with Sniffer data is not specified!", []);
    return;
elseif ( isempty(this.geda_param.snifrange) )
    this.make_report("dat", "WARNING: The range for Sniffer data is not specified!", []);
    return;
elseif ( isempty(this.geda_param.robid) )
    this.make_report("dat", "WARNING: The device name (robot) of Sniffer data is not specified!", []);
    return;
end

this.device = this.geda_param.robid;
this.outpath = "";
if ~isempty(this.geda_param.outpath)
    this.outpath = this.geda_param.outpath;
end

end


