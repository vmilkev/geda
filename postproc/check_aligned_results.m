
clear;


main("test/56593/param_56593_R6_2022_7.txt");
%%

p_file1 = "test/param_21834.txt"; % all good results
p_file1_2 = "test/param_21834_2.txt";
p_file1_3 = "test/param_21834_3.txt";
%--------------------------------------
p_file4 = "test/param_33489.txt"; % 10 sec spell exists
p_file4_2 = "test/param_33489_2.txt";
p_file4_3 = "test/param_33489_3.txt";
%--------------------------------------
p_file2 = "geda_param_105_v2.dat"; % 135 sec spell exists <= good improvement
p_file3 = "test/param_26089_2.txt"; % !!! check zeros in the output. Results are good.
%
p_file5 = "test/param_31563.txt"; % not enough memory
p_file6 = "test/param_57112.txt"; % good results. !!! check zeros in the output; 11-14 sec spells exist; NOTE, good reliability detection. !Denoising improves the detection.
p_file7 = "test/param_59551.txt"; % nice example of partly unreliable data!!!; 184 sec spell exists. Denoising improves the detection.
p_file8 = "test/param_63204.txt"; % all good results
p_file9 = "test/param_63625.txt"; % 1/3 of data not aligned; strange data, consider for fault detection case.
p_file10 = "test/param_110038.txt"; % 78-228 min spells exist <= no improvement. !!! Because the signalls are too short due to breacks.
p_file11 = "test/param_13144.txt"; % long running ...
p_file12 = "test/param_25184.txt"; % NOTE, good reliability detection
p_file13 = "test/31563_first/Handling_final/param_31563F_R2.txt";
p_file14 = "test/31563_second/Handling_final/param_31563S_R2.txt";
p_file15 = "test/CPH/Handling_test_1day/param_MO_1day_R1.txt";
p_file16 = "param_56614_R1.txt";
p_file17 = "test/param_84544_R7_2023_10.txt"; % partly reliable; moderately long running

% main(p_file6);

main(p_file1);
% main(p_file1_2);
% main(p_file1_3);

% main(p_file4);
% main(p_file4_2);
% main(p_file4_3);

%%
res = read_aligned("sniffer_R5_2022_7/aligned_18-Jun-2025_14-18-41_robot_5.geda"); % check results
%res = read_aligned("aligned_init_sniffer_4.txt"); % check results
%res = read_aligned("aligned_init_sniffer_2.txt"); % check results

%%

gases = size(res,2);

ch4 = {};
co2 = {};

get_ch4 = true;

figure(3);

gases = 5;

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

%% TEST SSA

t5 = duration(0,0,5); % 5 sec, knowing sampling interval
%tdif = (co2{1, 2}(2:end)-co2{1, 2}(1:end-1))-t5;
tdif = (ch4{1, 2}(2:end)-ch4{1, 2}(1:end-1))-t5;

testdata = table2array( res(:,4));
testdata = testdata(1:5000,1);

N = 200; % The number of time 'moments' in our toy series
N = 5000;

L = 70;
K = N - L + 1;

t = (1:N)';
F = testdata;

% making Hankel matrix
H = zeros(L,K);
for i1 = 1:L
    H(i1,1:K) = F( i1 : i1 + K-1 );
end

Hsqv = H*H';

[U,S,V] = svd(H);

% [coeff,score,latent,~,explained,mu] = pca(H);
% pca_95 = find(cumsum(explained)>50,1); % components explain more than 95% of all variability
% H95pca = coeff(:,1:pca_95)*coeff(:,1:pca_95)'*H;

noise_edge = 90; % procent of explained_eigenvalues
noise_edge80 = 80; % procent of explained_eigenvalues

explained_svd = 100*diag(S.^2)./sum(diag(S.^2));
svd_expl = find(cumsum(explained_svd)>noise_edge,1); % components explain more than 90% of all variability
svd_expl2 = find(cumsum(explained_svd)>noise_edge80,1); % components explain more than 80% of all variability

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
title("explained_svd'");
subplot(1,2,2);
plot(cumsum(explained_svd));
title("cumsum(explained_svd)'");

figure(6);
plot(t,F);
hold on;
for i = 1:svd_expl
    plot(t,RC(:,i));
end
hold off;
title("original (F) and components (RCi, 90 % variance) data");

figure(7);
plot(t,F, 'ro:', 'LineWidth', 2);
hold on;
y = sum(RC(:,1:svd_expl),2);
y2 = sum(RC(:,1:svd_expl2),2);
%plot(t,y, 'b');
plot(t,y2, 'k', 'LineWidth', 2);
hold off;
ttl = strcat("original (F) and reconstructed (RC, using ", num2str(svd_expl), " components) data");
title(ttl);

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
title("W");
