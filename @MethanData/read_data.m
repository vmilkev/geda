function read_data( this )
% Reads .xlsx & .csv data files

% we have to read Robot Data only once
if ( ~this.isLelyData )
    this.read_lely_data();
end

% here we can read mesurements data files multiple times
if ( ~isempty(this.geda_param) )
    this.read_mesr_data();
end

end
