function r_name = get_mesr_robot( this )
% Public inteface to the class data member this.mesr_param.robid
    r_name = this.mesr_param.robid;
    if isempty( r_name )
        disp("get_mesr_robot(): There is no known robot ID!");
    end
end
