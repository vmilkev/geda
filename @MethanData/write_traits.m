function write_traits( this, summary_table, ams_fnames, snf_fnames, used_sig_length )

trait_data_file_1 = strcat("trait_v1_",...
    this.date_now,...
    "_robot_",...
    num2str(this.device),... % robot (device) id
    ".geda");

if this.outpath ~= ""
    if ~exist(this.outpath, 'dir')
        mkdir(this.outpath);
    end
    trait_data_file_1 = strcat(this.outpath,trait_data_file_1);
end

if exist(trait_data_file_1, 'file')
    this.make_report("dat", "WARNING: There is the trait file of the same name already exist in the output folder!", []);
    this.make_report("dat", "         The existing data will be overwriten by the new data!", []);
    delete(trait_data_file_1);
end

this.make_report("dat", strcat("Calculating traits and writing into the file: ", trait_data_file_1), []);

varNames = {'time','id'};

ams_ts = this.get_from_binary( ams_fnames{1} ); % loading AMS TS with 1 sec frequency; assuming only one AMS TS

Fts = this.sampl_res_const;
Fact = this.deltaT_original;
i_signal = 0;
bkg_outl_const = -0.1; % background outlier
obs_outl_const = 1.5; % observations outlier in statistic: data./(3*std(data)) > obs_outl_const -> is outlier if TRUE;
wlen_const = 10; % SSA window, sec
explvar_const = 90;

use_cols = true(1,1); % preallocating logical array indicating gas cols used for making traits
rlb_list = find( summary_table(:,7) == 1 ); % only reliable data

all_traits = cell(numel(rlb_list), 1);

for i = 1:size(rlb_list,1) % loop over the number of processed signals

    if summary_table( rlb_list(i), 6 ) ~= i_signal % only when we change the dataset
        i_signal = summary_table( rlb_list(i), 6 ); % the num id of the dataset stored in the binary file
        snf_ts = this.get_from_binary(snf_fnames{i_signal,1}); % load Sniffer TS with 1 sec sampling frequency
    end

    first_ams = ( summary_table( rlb_list(i), 2 ) - 1 ) * Fts + 1; % first index of the data range in ams_ts
    last_ams = ( summary_table( rlb_list(i), 3 ) - 1 ) * Fts + 1; % last index of the data range in ams_ts
    first_snf = ( summary_table( rlb_list(i), 4 ) - 1 ) * Fts + 1; % first index of the data range in snf_ts
    last_snf = ( summary_table( rlb_list(i), 5 ) - 1 ) * Fts + 1; % last index of the data range in snf_ts

    % create resampling indeces
    l = 1;
    for j = first_ams:Fact:last_ams
        samples_ams(l,1) = uint32(j);
        l = l + 1;
    end
    l = 1;
    for j = first_snf:Fact:last_snf
        samples_snf(l,1) =  uint32(j);
        l = l + 1;
    end

    t = datetime( ams_ts( samples_ams,1 ),'ConvertFrom','datenum','Format','dd-MMM-yyy HH:mm:ss' );
    id = uint64( ams_ts( samples_ams,3 ) );
    gas = snf_ts(samples_snf,2:end);

    %------------------------------------
    if i == 1 % do this only once
        % select meaningfull cols (indicating that something's changing over time)
        right_cols = zeros(size(gas,2), 1);
        for j = 1:size(gas,2)
            if (std(gas(:,j)) == 0) || any(isnan(gas(:,j)))
                continue;
            else
                right_cols(j,1) = j;
            end
        end
        use_cols = logical(right_cols);
    end

    gas = gas(:,use_cols);

    [t,id,gas] = filter_baddata(t,id,gas);

    mask_bkg = id == 0; % mask the background dataset

    gas_bkg = gas(mask_bkg,:);

    bkg_val = zeros(size(gas,2),1); % background concentrations
    for j = 1:size(gas,2)
        bkg_val(j,1) = mode(gas_bkg(gas_bkg(:,j)>bkg_outl_const,j)); % the most frequent values above the outlier level at no-measurement period
        gas(:,j) = gas(:,j) - bkg_val(j,1); % correct for background
        f = gas(:,j) < 0.0; % filtering records below the background
        gas(f,j) = 0.0; % this is supposed to be the corrected background
    end

    clear gas_bkg
    
    Hgas = gas(~mask_bkg,:);
    Hid = id(~mask_bkg,:);
    Htime = t(~mask_bkg,:);

    clear samples_ams samples_snf gas t id mask_bkg

    for l = 1:size(Hgas,2)
        std_val = 3.0 * std( Hgas(:,l) ); % outlier upper bound
        f = Hgas(:,l) ./ std_val > obs_outl_const; % find outliers
        Hgas(f,l) = 0.66 * std_val; % correct outliers
        [ H1,~,~,~ ] = nstd_ssa( Hgas(:,l), wlen_const, [], explvar_const );
        f2 = H1 < 0.0; % filtering records below the background
        H1(f2) = 0.0; % this is supposed to be the corrected background
        Hgas(:,l) = H1;
    end

    clear H1 f
    
    u = unique(Hid);
    
    n_traits = numel(u);
    jstart = 1;
    if u(1) == 0
        jstart = 2;
        n_traits = numel(u) - 1;
        this.make_report("datint", "WARNING: The 0 ID present at the measurements interval. Reliable signal at pos.: ", i);
    end

    traits = zeros(n_traits,size(Hgas,2));
    tr_ids = zeros(n_traits,1);
    tr_time = strings(n_traits,1);
    background = zeros(n_traits,size(Hgas,2));
    
    j2 = 0;
    for j = jstart:numel(u)
        animal = u(j);
        fnd_id = Hid == animal;
        anim_gas = Hgas(fnd_id,:);
        anim_time = Htime(fnd_id,:);

        t_duration = numel(anim_gas(:,1))*Fact/60.0; % duration in min

        if t_duration < 2.0
            continue;
        end

        j2 = j2 + 1;

        for l = 1:size(Hgas,2)
            traits(j2,l) = trapz( Fact, anim_gas(:,l) )/t_duration;
            background(j2,l) = bkg_val(l,1);
        end

        tr_ids(j2,1) = animal;
        tr_time(j2,1) = anim_time(1,1);
    end

    clear Htime Hid Hgas anim_gas anim_time fnd_id

    T = table( tr_time(1:j2,1), tr_ids(1:j2,1), 'VariableNames', varNames );
    for i_gas = 1:2
        i_name1 = strcat("gas_", num2str(i_gas));
        i_name2 = strcat("bkg_", num2str(i_gas));
        T.(i_name1) = traits(1:j2,i_gas);
        T.(i_name2) = background(1:j2,i_gas);
    end
    
    all_traits{i,1} = T;

    clear T

end

writetable(vertcat(all_traits{:}),trait_data_file_1,'WriteMode','Append',...
    'WriteVariableNames',true,'Delimiter',';',...
    'QuoteStrings',true, 'FileType', 'text');

[ch4, co2] = fit_model(all_traits,used_sig_length);

clear all_traits

this.make_report("dat", " ", []);
this.make_report("dat", "Summary of LMM analysis of the processed data:", []);
this.make_report("dat", strcat("        CH4 model: ", ch4{1,1}), []);
this.make_report("dat", strcat("        num. observations per ID:  ", num2str(ch4{4,1}) ), []);
this.make_report("dat", strcat("        fit quality, adjusted r^2: ", num2str(ch4{3,1}) ), []);
this.make_report("dat", strcat("        var_id/(var_id+var_err):   ", num2str(ch4{2,1}) ), []);
this.make_report("dat", " ", []);
this.make_report("dat", strcat("        CO2 model: ", co2{1,1}), []);
this.make_report("dat", strcat("        num. observations per ID:  ", num2str(co2{4,1}) ), []);
this.make_report("dat", strcat("        fit quality, adjusted r^2: ", num2str(co2{3,1}) ), []);
this.make_report("dat", strcat("        var_id/(var_id+var_err):   ", num2str(co2{2,1}) ), []);
this.make_report("dat", " ", []);

clear ch4 co2

% ---------------------------------------------------------------
% Nested Functions
% ---------------------------------------------------------------
    function [flt_t, flt_id, flt_gas] = filter_baddata(time, ids, gases)
        mask = true(size(gases,1),1);
        for igas = 1:size(gases,2) % loop over gas types
            mask_nan = []; % mask NAN records
            if sum( isnan(gases(:,igas)) ) > 0
                mask_nan = ~isnan(gases(:,igas));
            end
            mask_inf = []; % mask INF records
            if sum( isinf(gases(:,igas)) ) > 0
                mask_inf = ~isinf(gases(:,igas));
            end
            if ~isempty(mask_nan)
                mask = mask & mask_nan;
            end
            if ~isempty(mask_inf)
                mask = mask & mask_inf;
            end
        end
        flt_gas = gases(mask,:);
        flt_t = time(mask,:);
        flt_id = ids(mask,:);
    end

    function [ H,l,nl,F ] = nstd_ssa( X, L, Q, E )
        % X := signal to analyse
        % L := size of window for Hankel matrix
        % Q := number of eigen components to use in reconstruction, n <= L;
        %      if Q == [], Q will be defined based on argument E
        % E := procent of variance explained by eigenvalues;
        %      this value defines number of eigenvalues (Q)
        %      used for reconstruction; this option is active if Q == []
        %
        % H := reconstructed signal
        % l := vector of (all) eigenvalues, dim(l) = dim(L)
        % nl := number of eigenvalues responsible for E % variation
        % F := first nl number of eigenvectors

        N = size(X,1);

        Y=zeros(N-L+1,L);
        for m=1:L
            Y(:,m) = X((1:N-L+1)+m-1);
        end

        Cemb=Y'*Y / (N-L+1);
        C=Cemb;

        [V,LAMBDA] = eig(C);
        LAMBDA = diag(LAMBDA);               % extract the diagonal elements
        [LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
        V = V(:,ind);                    % and eigenvectors

        PC = Y*V;

        RC=zeros(N,L);
        for m=1:L
            buf=PC(:,m)*V(:,m)'; % invert projection
            buf=buf(end:-1:1,:); % revert rows, in order to work with subdiagonals
            for n=1:N % anti-diagonal averaging
                RC(n,m)=mean( diag(buf,-(N-L+1)+n) ); % Get the elements on the (N-L+1)+n subdiagonal (k=-(N-L+1)+n) of A.
            end
        end

        if ~exist('E','var') || isempty(E)
            if isempty(Q)
                disp('Error: At least one of these two parameters should be supplied: Q or E!');
                return
            else
                H = sum(RC(:,1:Q),2);
                nl = Q;
            end
        else
            l_sum = cumsum(LAMBDA)*100./sum(LAMBDA);
            n_E = find( l_sum <= E );
            if isempty(n_E)
                nl = 1;
            else
                nl = max(n_E) + 1;
            end
            H = sum(RC(:,1:nl),2);
        end

        l = LAMBDA;
        F = V;
    end

    function [res_ch4, res_co2] = fit_model(trait_list, used_sig_len)           
        
        tr = vertcat(trait_list{:});
        clear trait_list
        
        u_tr = unique(tr.id);
        obs_per_id = size(tr,1)/numel(u_tr);
        
        tn = convertTo(datetime(tr.time), 'posixtime'); 
        tn = round((tn - min(tn) + 1)./60);
        tc = zeros( size(tn) );
        h = used_sig_len; % used sig. length
        if used_sig_len > 23
            h = 1;
        end
        period = h*60; % observation period in 24 hours (n_perionds = (24/h)*60), min
        day = 24*60;
        for i_tr = 1:size(tn,1)
            rel_time = floor((tn(i_tr,1) - floor(tn(i_tr,1)/day)*day)/period);
            tc(i_tr,1) = rel_time + 1;
        end
        
        tday = datetime(tr.time);
        tday.Format = 'dd-MMM-yyyy';
        tday = categorical(string(tday));
        
        tr.("period") = tc;
        tr.("day") = tday;
                
        model_ch4 = 'gas_1 ~ 1 + period + day + bkg_1 + (1|id)';    
        model_co2 = 'gas_2 ~ 1 + period + day + bkg_2 + (1|id)';
        
        lme = fitlme(tr, model_ch4);
        [psi,mse] = covarianceParameters(lme);
        h_ch4 = psi{1}/(psi{1}+mse); % repeatability/heritability
        r_sqv_ch4 = lme.Rsquared.Adjusted;
        res_ch4{1,1} = model_ch4;
        res_ch4{2,1} = h_ch4;
        res_ch4{3,1} = r_sqv_ch4;
        res_ch4{4,1} = obs_per_id;
        
        lme2 = fitlme(tr, model_co2);
        [psi,mse] = covarianceParameters(lme2);
        h_co2 = psi{1}/(psi{1}+mse);
        r_sqv_co2 = lme2.Rsquared.Adjusted;
        res_co2{1,1} = model_co2;
        res_co2{2,1} = h_co2;
        res_co2{3,1} = r_sqv_co2;
        res_co2{4,1} = obs_per_id;       
    end
% ---------------------------------------------------------------

end