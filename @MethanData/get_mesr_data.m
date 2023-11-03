function data = get_mesr_data( this, type )

% PUBLIK INTERFACE TO THE CLASS MEMBER VARIABLES
% RETURNS A MESUREMENTS DATA

if ~this.isMesrData
    this.read_data();
end

if isempty(this.MesrMap)
    disp("get_mesr_data(): The this.MesrMap data is empty.");
    return;
end

data = cell( length(this.MesrMap),1 );

if nargin < 2 % ............................ Return the NUMERIC data as the defoult option
    for i = 1:length(this.MesrMap)
        data{i,1} = this.MesrNumMap(i);
    end
    return;
end

if ( nargin == 2 )
    if strcmp("numeric", type) % ........... Return the NUMERIC data
        for i = 1:length(this.MesrMap)
            data{i,1} = this.MesrNumMap(i);
        end
    elseif strcmp("table", type) % ......... Return TABLE data
        for i = 1:length(this.MesrMap)
            data{i,1} = this.MesrMap(i);
        end
    else
        disp("get_mesr_data(): Wrong input argument.");
    end
end

end