
clear;

p_file1 = "test/param_21834.txt";
p_file2 = "geda_param_105_v2.dat"; %                               135 sec spell exists <= good improvement
p_file3 = "test/param_26089.txt"; % !!! check zeros in the output. Results are good.
p_file4 = "test/param_33489.txt"; %                                10 sec spell exists
p_file5 = "test/param_31563.txt"; % not enough memory
p_file6 = "test/param_57112.txt"; % good results. !!! check zeros in the output; 11-14 sec spells exist; NOTE, good reliability detection. !Denoising improves the detection.
p_file7 = "test/param_59551.txt"; % nice example of partly unreliable data;      184 sec spell exists. Denoising improves the detection.
p_file8 = "test/param_63204.txt"; % good results
p_file9 = "test/param_63625.txt"; % 1/3 of data not aligned; strange data, consider for fault detection case.
p_file10 = "test/param_110038.txt"; %                              78-228 min spells exist <= no improvement. !!! Because the signalls are too short due to breacks.
p_file11 = "test/param_13144.txt"; % long running ...
p_file12 = "test/param_25184.txt"; % NOTE, good reliability detection
p_file13 = "test/31563_first/Handling_final/param_31563F_R2.txt";
p_file14 = "test/31563_second/Handling_final/param_31563S_R2.txt";
p_file15 = "test/CPH/Handling_test_1day/param_MO_1day_R1.txt";
p_file16 = "param_56614_R1.txt";

tic;
main(p_file1); % run alignment
toc

%%
[res] = read_aligned("aligned_adj_sniffer_103.txt"); % check results

%%

gases = size(res,2);

ch4 = {};
co2 = {};

get_ch4 = true;

figure(3);

for i_plot = 3:gases-1

    i_var = res.Properties.VariableNames{i_plot};

    which_gas = table2array( res(:,i_plot) );
    
    id = table2array(res(:,2));
    u = unique(id); % sorted order, zero comes first
    
    min_val = min( which_gas(which_gas > 0.0) );
    if (isempty(min_val))
        break;
    end

    id( id == 0 ) = min_val;
    
    mean_val = mean( which_gas(which_gas > 0.0) );
    std_val = std( which_gas(which_gas > 0.0) );

    all_time = table2array( res(:,1) );
    
    % Scaling animals' IDs
    % starting from 2 because the first ID is zero
    for i = 2:size(u,1)
        
        j = find(id == u(i,1));
        %id(j) = 1+i*0.005;
        id(j) = mean_val+rand*sqrt(std_val);
        
        if (get_ch4)
            ch4{i-1,1} = u(i,1);
            co2{i-1,1} = u(i,1);
            records = table2array( res(j,3) );
            time = table2array( res(j,1) );
            mask = records > 0;
            ch4{i-1,3} = records(mask);
            ch4{i-1,2} = time(mask);
            
            records = table2array( res(j,4) );
            co2{i-1,3} = records;
            co2{i-1,2} = time;
        end

    end

    if (get_ch4)
        get_ch4 = false;
    end
        
    subplot(gases-3,1,i_plot-2);
    plot( all_time, which_gas, 'o', all_time, id, 'o' );
    xlabel( 'time' );
    ylabel(strcat('gas no.',num2str(i_plot-2)));

end

%% TEST PCA (SSA)
t5 = duration(0,0,5); % 5 sec, knowing sampling interval
%tdif = (co2{1, 2}(2:end)-co2{1, 2}(1:end-1))-t5;
tdif = (ch4{1, 2}(2:end)-ch4{1, 2}(1:end-1))-t5;
stop_point = find(tdif > duration(0,0,0));

which_set = 1;

p1 = stop_point(which_set)+1;
p2 = stop_point(which_set+1)-1;
%testdata = co2{1, 3}(p1:p2);
%testdata = ch4{1, 3}(p1:p2);
%testdata = co2{1, 3}(:);
testdata = table2array( res(:,4));% ids = table2array( res(:,2) ); mask = ids > 0; testdata = testdata(mask,1); testdata = testdata(1:5000,1);
testdata = testdata(1:5000,1);

%%
%--------------------------------------------------------------------
N = 200; % The number of time 'moments' in our toy series
N = 5000;

L = 70;
K = N - L + 1;

t = (1:N)';
F = testdata;

% trend = 0.001 * (t - 100).^2;
% p1 = 20;
% p2 = 30;
% periodic1 = 2 * sin(2*pi*t/p1);
% periodic2 = 0.75 * sin(2*pi*t/p2);
% 
% rng(123) % So we generate the same noisy time series every time.
% noise = 2 * (rand(N,1) - 0.5);
% F = trend(:,1) + periodic1(:,1) + periodic2(:,1) + noise(:,1);
% 
% figure(1);
% clf;
% hold on;
% plot(t,F);
% plot(t,trend);
% plot(t,periodic1);
% plot(t,periodic2);
% plot(t,noise);
% hold off;
%--------------------------------------------------------------------


% making Hankel matrix
H = zeros(L,K);
for i1 = 1:L
    H(i1,1:K) = F( i1 : i1 + K-1 );

end

figure(2);
imagesc(H);


Hsqv = H*H';

figure(3);
imagesc(Hsqv);

[U,S,V] = svd(H);

% [coeff,score,latent,~,explained,mu] = pca(H);
% pca_95 = find(cumsum(explained)>50,1); % components explain more than 95% of all variability
% H95pca = coeff(:,1:pca_95)*coeff(:,1:pca_95)'*H;

noise_edge = 90; % procent of explained_eigenvalues

explained_svd = 100*diag(S.^2)./sum(diag(S.^2));
svd_expl = find(cumsum(explained_svd)>noise_edge,1); % components explain more than 95% of all variability

% Reconstructing
RC=zeros(N,L);
for m=1:L
  buf = S(m,m)*U(:,m)*V(:,m)';
  buf=buf(end:-1:1,:); % revert rows, in order to work with subdiagonals
  for n=1:N % anti-diagonal averaging
    ind = -(L-1)+n-1;
    RC(n,m)=mean( diag(buf,ind) ); % Get the elements on the (N-L+1)+n subdiagonal (k=-(N-L+1)+n) of A.
  end
end

figure(5);
subplot(1,2,1);
plot(explained_svd);
subplot(1,2,2);
plot(cumsum(explained_svd));

figure(6);
plot(t,F);
hold on;
for i = 1:12
    plot(t,RC(:,i));
end
hold off;

figure(7);
plot(t,F);
hold on;
plot(t,sum(RC(:,1:4),2));
hold off;

w = zeros(N,1);
for i = 1:N
    if (i >= 1 && i <= L)
        w(i,1) = i + 1;
    elseif (i >= L+1 && i <= K)
            w(i,1) = L;
    elseif (i >= K+1 && i <= N)
            w(i,1) = N - i;
    end
end

W = zeros(L);
for i = 1:L
    for j = 1:L
        W(i,j) = ( RC(:,i)'*(RC(:,j).*w) )/( sqrt((RC(:,i)'*(RC(:,i).*w))) * sqrt((RC(:,j)'*(RC(:,j).*w))) );
    end
end

figure(8);
imagesc(W);


%% DYNAMICAL NOISE

which_gas = 4;
range = [1:70];
ids = table2array( res(range,2) );
records = table2array( res(range,which_gas) );
time = table2array( res(range,1) );
mask = ids > 0;
all_gas = records(mask,1);
time_gas = time(mask,1);
% all_gas = records;
% time_gas = time;


Lmk = {};
N = size(all_gas,1);
all_k = 500;

for i_k = 1:all_k
    for i_m = 1:i_k
        %L = zeros(1,floor( (N-i_m)/i_k )+1);
        for i = 0:floor( (N-i_m)/i_k )
            L(i+1) = all_gas(i_m+i*i_k,1);
        end
        Lmk{i_k,1}(i_m) = (sum(abs(L(2:end)-L(1:end-1)))*( (N-1) / (i_k*floor( (N-i_m)/i_k )) ))/i_k;
    end
end


for i = 1:all_k
    Lmean(i,1) = mean(Lmk{i,1}(:));
end 


k2 = (1./[1:all_k]);

k3 = log(k2);
L3 = log(Lmean);

figure;
plot((k2),log(Lmean),'o-');

%% CREATE TIMETABLE FOR THE SPECIFIC ID

which_data = 1;
i_data = datenum(ch4{which_data,2});
i_data = (i_data - min(i_data)) *24*60*60; % transform time to sec starting from zero
s_data = timetable(seconds(i_data(403:477,1)),ch4{which_data,3}(403:477,1));
%s_data = timetable(seconds(0:920)',ch4{which_data,3});
t = seconds(uint64(i_data));
v = ch4{which_data,3};
% figure;
% nplots = 5;
% first = 100;
% vMax1 = 0;
% for i = 1:nplots
%     subplot(nplots,1,i);
%     j = find(res.id == u(first+1,1));
%     first = first + 1;
%     plot(res.t(j), res.ch4(j),'.');
%     vMax2 = max(res.ch4(j));
%     if (vMax2 > vMax1)
%         vMax1 = vMax2;
%     end
%     ylim([0 vMax1]);
% end

% ind_f = find( abs(q(:,2)) > 200);
% ind_t = find( abs(q(:,2)) < 200);
% q_f = q( ind_f,: );
% q_t = q( ind_t,: );
% figure, hold on, plot(q_t(:,5), q_t(:,7), 'bo'), plot(q_f(:,5), q_f(:,7), 'ro')

%%

function [res] = read_aligned(filename, dataLines)

%     % Input handling
%     
%     % If dataLines is not specified, define defaults
%     if nargin < 2
%         dataLines = [2, Inf];
%     end
%     
%     %% Set up the Import Options and import the data
%     opts = delimitedTextImportOptions("NumVariables", 4);
%     
%     % Specify range and delimiter
%     opts.DataLines = dataLines;
%     opts.Delimiter = ";";
%     
%     % Specify column names and types
%     opts.VariableNames = ["t", "id", "co2", "ch4"];
%     opts.VariableTypes = ["datetime", "double", "double", "double"];
%     
%     % Specify file level properties
%     opts.ExtraColumnsRule = "ignore";
%     opts.EmptyLineRule = "read";
%     
%     % Specify variable properties
%     opts = setvaropts(opts, "t", "InputFormat", "dd-MMM-yyyy HH:mm:ss");
%     
%     % Import the data
%     res = readtable(filename, opts);

    if nargin < 2
        dataLines = [2, Inf];
    end

    opts = detectImportOptions(filename);

    % Specify range
    opts.DataLines = dataLines;
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    opts = setvaropts(opts, opts.VariableNames{1}, "InputFormat", "dd-MMM-yyyy HH:mm:ss");
    
    % Import the data
    res = readtable(filename, opts);

end