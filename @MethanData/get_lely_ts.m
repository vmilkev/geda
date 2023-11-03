function data = get_lely_ts( this, robots, delta_t )

    % PUBLIK INTERFACE TO THE CLASS MEMBER VARIABLES
    % RETURNS A ROBOT TIME SERIES DATA
    
    if ~this.isLelyData
        this.read_data();
    end
    
    if isempty(this.Robots)
        disp("get_lely_ts(): The ROBOT data is empty.");
        return;
    end
    
    data = cell( numel(robots),1 );
    
    for i = 1:numel( robots )
        data{i,1} = this.create_ts2( this.LelyNumMap(robots(i)), 1.0 ); % extend to 1 sec interval
        if ( delta_t > 1.0 )
            data{i,1} = this.shrink( data{i,1}, delta_t ); % shrink to requested deltaT
        end
    end

end