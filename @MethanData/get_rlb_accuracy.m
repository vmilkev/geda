function [accuracy, corrected_rlb] = get_rlb_accuracy( this, rlb_estimates )

    limit_const = 47000; % 12 hours + 10%, in [sec]
    f0 = find ( abs(rlb_estimates(:,1)) > limit_const ); % supposed to be the reference false => 0
    f1 = find ( abs(rlb_estimates(:,1)) <= limit_const ); % supposed to be the reference true => 1
    
    f0_wrong = find( rlb_estimates(f0,2) == 1 ); % incorrectly detected true
    f1_wrong = find( rlb_estimates(f1,2) == 0 ); % incorrectly detected false
    
    accuracy = 1.0 - ( numel(f0_wrong) + numel(f1_wrong) )/size(rlb_estimates,1);
    
    corrected_rlb = zeros( size(rlb_estimates,1),1 );
    
    if accuracy < 1.0
        corrected_rlb(f0,1) = 0;
        corrected_rlb(f1,1) = 1;
    else
        corrected_rlb(:,1) = rlb_estimates(:,2);
    end

end