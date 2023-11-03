function [data, robots] = get_lely_data( this, type )

% PUBLIK INTERFACE TO THE CLASS MEMBER VARIABLES
% RETURNS A ROBOT DATA

if ~this.isLelyData
    this.read_data();
end

if isempty(this.Robots)
    disp("get_lely_data(): The ROBOT data is empty!");
    return;
end

data = cell( numel(this.Robots),1 );

if nargin < 2 % ........................................... Return the NUMERIC data as the defoult option
    for i = 1:numel( this.Robots )
        data{i,1} = this.LelyNumMap( this.Robots(i) );
    end
    robots = this.Robots;
    return;
end

if ( nargin == 2 )
    if strcmp("numeric", type) % .......................... Return the NUMERIC data
        for i = 1:numel( this.Robots )
            data{i,1} = this.LelyNumMap( this.Robots(i) );
        end
    elseif strcmp("table", type) % ........................ Return TABLE data
        for i = 1:numel( this.Robots )
            data{i,1} = this.LelyMap( this.Robots(i) );
        end
    else
        disp("get_lely_data(): Wrong input argument.");
    end
    robots = this.Robots;
end


end