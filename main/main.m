function main( prmfile )

    % prmfile - the main parameter file where:
    %           1) names & paths to AMS param file and Sniffer param files are specified;
    %           2) parameters defined, such as the sampling interval and skews estimates. 
    
    warning('off','all');
    
    if nargin < 1
        disp("There is no parameter file!");
        return;
    end
        
    alignTS( prmfile );

end