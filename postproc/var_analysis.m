clear;
close all;

%make_phenotypes("aligned_init_sniffer_4.txt",4, "ids_list_sniffer_4.txt");
%make_phenotypes("aligned_init_sniffer_2.txt", 4, "ids_list_sniffer_2.txt");
make_phenotypes("aligned_rlb_sniffer_4.txt",4);
%make_phenotypes("aligned_init_sniffer_1.txt",4, "ids_list_sniffer_1.txt");
%make_phenotypes("aligned_init_sniffer_101.txt",4);

%%
device = ["1" "2" "3" "4"];
%device = ["1" "2"];
%device = ["3" "4"];
gas = "1";
type = "rlb";
%type = "init";

files = {};
template = ["traits_gas_" "_aligned_" "_sniffer_" ".txt"];
for i = 1:size(device,2)
    files{i} = strcat(template(1), gas, template(2), type, template(3),device(i),template(4));
end

tn = [];
idn = [];
obs = [];
dev = [];
u_ids = {};

for i = 1:size(device,2)
    data = import_traits(files{i});
    tn_ = convertTo(data.time, 'posixtime');   
    idn_ = data.id;

    u_ids{i} = unique(idn_);
    
    dev_ = zeros(size(tn_));
    dev_(:,1) = str2num(device(i));
    
    disp(["Unique ids:" num2str(numel(unique(idn_)))]);
    
    % observations
    obs1 = data.trait_1;
    obs2 = data.trait_2;
    obs3 = data.trait_3;
    obs4 = data.trait_4;
    obs_ = [obs1, obs2, obs3, obs4];

    % -----------------------------
    % remove outlier    
    % obs(66, :) = 0;
    % mask = obs(:, 1) > 0;
    % obs = obs(mask,:);
    % tn = tn(mask,:);
    % tc = tc(mask,:);
    % idn = idn(mask,:);
    % -----------------------------

    tn = [tn; tn_];
    idn = [idn; idn_];
    obs = [obs; obs_];
    dev = [dev; dev_];
end

disp(["Overall (combined) unique ids:" num2str(numel(unique(idn)))]);

% find ids common to all datasets
% u = unique(idn);
% common_ids = zeros(size(u));
% for i = 1:size(u,1)
%     i_id = u(i);
%     is_comm = ismember(i_id, u_ids{1});
%     for j = 2:2 %size(u_ids,2)
%         is_comm2 = ismember(i_id, u_ids{j});
%         is_comm = is_comm & is_comm2;
%     end
%     if is_comm
%         common_ids(i,1) = i_id;
%     end
% end
% mask = find(common_ids ~= 0);
% common_ids = common_ids(mask);
% mask = idn == common_ids(1,1);
% for i = 2:size(common_ids,1)
%     mask2 = idn == common_ids(i,1);
%     mask = mask + mask2;
% end

tn = round((tn - min(tn) + 1)./60);

tc = zeros( size(tn) );
h = 4;
period = h*60; % observation period in 24 hours (n_perionds = (24/h)*60), min
day = 24*60;
for i = 1:size(tn,1)
    rel_time = floor((tn(i,1) - floor(tn(i,1)/day)*day)/period);
    tc(i,1) = rel_time + 1;
end

tb = table(tn, tc, dev, idn, obs(:,1), obs(:,2), obs(:,3), obs(:,4), 'VariableNames',{'time', 'period', 'device', 'id', 'obs1', 'obs2', 'obs3', 'obs4'});

lme1 = fitlme(tb, 'obs1 ~ 1 + period + device + (1|id)');
lme2 = fitlme(tb, 'obs2 ~ 1 + period + device + (1|id)');
lme3 = fitlme(tb, 'obs3 ~ 1 + period + device + (1|id)');
lme4 = fitlme(tb, 'obs4 ~ 1 + period + device + (1|id)');

lme{1} = lme1;
lme{2} = lme2;
lme{3} = lme3;
lme{4} = lme4;
    
figure(1)
set(gcf,'Color','white')

for i = 1:4
    [psi,mse] = covarianceParameters(lme{i});
    if i == 1
        disp(["sigmma", "err", "h", "R^2"]);
    end
    disp([psi{1}, mse, psi{1}/(psi{1}+mse), lme{i}.Rsquared.Adjusted]);
    
    [B,Bnames,stats] = randomEffects(lme{i});
    % ---------------------
    iobs = obs(:,i);
    yhat = predict(lme{i},tb);
    subplot(1,4,i);
    %hold on
    plot(iobs, yhat, '*');
    xlim([min(iobs) max(iobs)]);
    ylim([min(iobs) max(iobs)]);
    box on
    axis square
    % ---------------------
end

u = unique(idn);
ids_var = zeros(size(u,1), 4);
ids_std = zeros(size(u,1), 4);
ids_mean = zeros(size(u,1), 4);
for i = 1:size(u,1)
    mask = idn == u(i);
    for j = 1:size(obs,2)        
        ids_var(i,j) = var( obs(mask, j) );
        ids_std(i,j) = std( obs(mask, j) )/mean(obs(mask, j));
        ids_mean(i,j) = mean( obs(mask, j) );
    end
end

all_var = zeros(1, size(obs, 2));
all_std = zeros(1, size(obs, 2));
across_var = zeros(1, size(obs, 2));
across_std = zeros(1, size(obs, 2));
ratio_std = zeros(1, size(obs, 2));
for i = 1:size(obs,2)
    all_var(1,i) = var(obs(:,i));
    all_std(1,i) = std(obs(:,i))/mean(obs(:,i));
    across_var(1,i) = var(ids_mean(:,i));
    across_std(1,i) = std(ids_mean(:,i))/mean(ids_mean(:,i));
    mask = ids_std(:,i) <= across_std(1,i);
    ratio_std(1,i) = numel(ids_std(mask,i))/numel(ids_std(:,i));
end

% figure(2);
% x = linspace(1, numel(u));
% y = zeros(size(x));
% y2 = zeros(size(x));
% for i = 1:size(obs, 2)
%     subplot(1,size(obs, 2),i);
%     hold on
%     plot(ids_var(:,i), '*');
%     y(:) = all_var(1,i);
%     y2(:) = across_var(1,i);
%     line(x,y,'Color','red','LineStyle','-', 'LineWidth', 2);
%     line(x,y2,'Color','black','LineStyle','-', 'LineWidth', 2);
%     xlim([0 max(x)]);
% end

figure(3);
set(gcf,'Color','white')
x = linspace(1, numel(u));
y = zeros(size(x));
y2 = zeros(size(x));
for i = 1:size(obs, 2)
    subplot(2,size(obs, 2),i);
    hold on
    plot(ids_std(:,i), '*');
    y(:) = all_std(1,i);
    y2(:) = across_std(1,i);
    line(x,y2,'Color','red','LineStyle','-', 'LineWidth', 2);
    line(x,y,'Color','black','LineStyle','-.', 'LineWidth', 1);
    ylim([0 max(max(ids_std))]);
    xlim([0 max(x)+1]);
    xlabel("ID's observations");
    ylabel("relative std");
    title(["variances below threshold: ", num2str(ratio_std(1,i))]);
    box on
end
for i = size(obs, 2)+1:size(obs, 2) * 2
    subplot(2,size(obs, 2),i);
    hold on
    histogram(ids_std(:,i-size(obs, 2)), 'Normalization','probability');
    ylabel("probability");
    xlabel("relative std");
    box on
end

%% FUNCTIONS

function d = import_traits(filename, dataLines)

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["time", "id", "trait_1", "trait_2", "trait_3", "trait_4"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "time", "InputFormat", "dd-MMM-yyyy HH:mm:ss");

% Import the data
d = readtable(filename, opts);

end