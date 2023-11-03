function [ sampT, siglen, is_ssa, ssa1, ssa2 ] = set_geda_param( this, in_param )

% Get necessary properties of LelyData file and set lely_param structure
% This is overloaded method. The IN argument can be
% 1) string/char array (file name) or 2) structure of lely_param type.

if ( isstring(in_param) || ischar(in_param) ) % check if IN argument is a file name
    this.geda_param = this.read_fparam( in_param );
elseif ( isstruct(in_param) ) % otherwise, check if it is a structure
    this.geda_param = in_param;
else
    disp("There is some problem with supplied parameter's TYPE/VALUE.");
    return;
end

sampT = this.geda_param.deltaT;
siglen = this.geda_param.siglength;
ssa1 = this.geda_param.wndlength;
ssa2 = this.geda_param.eigenvalues;
is_ssa = this.geda_param.denoise;

end


