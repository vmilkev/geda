clear;
close all;

make_phenotypes("aligned_init_sniffer_4.txt",4, "ids_list_sniffer_4.txt");
%make_phenotypes("aligned_init_sniffer_1.txt",4, "ids_list_sniffer_1.txt");
%make_phenotypes("aligned_init_sniffer_101.txt",4);

%%

data = import_traits("traits_gas_1_aligned_init_sniffer_4.txt");

tn = convertTo(data.time, 'posixtime');
tn = round((tn - min(tn) + 1)./60);

tc = zeros( size(tn) );
h = 4;
period = h*60; % observation period in 24 hours (n_perionds = (24/h)*60), min
day = 24*60;
for i = 1:size(tn,1)
    rel_time = floor((tn(i,1) - floor(tn(i,1)/day)*day)/period);
    tc(i,1) = rel_time + 1;
end

% IDs
idn = data.id;

u = unique(idn);
disp(["Unique ids:" num2str(numel(u))]);

% observations
obs1 = data.trait_1;
obs2 = data.trait_2;
obs3 = data.trait_3;
obs4 = data.trait_4;
obs = [obs1, obs2, obs3, obs4];
% -----------------------------
% remove outlier

% obs(66, :) = 0;
% mask = obs(:, 1) > 0;
% obs = obs(mask,:);
% tn = tn(mask,:);
% tc = tc(mask,:);
% idn = idn(mask,:);
% -----------------------------

tb = table(tn, tc, idn, obs(:,1), obs(:,2), obs(:,3), obs(:,4), 'VariableNames',{'time', 'period', 'id', 'obs1', 'obs2', 'obs3', 'obs4'});

lme1 = fitlme(tb, 'obs1 ~ 1 + period + (1|id)');
lme2 = fitlme(tb, 'obs2 ~ 1 + period + (1|id)');
lme3 = fitlme(tb, 'obs3 ~ 1 + period + (1|id)');
lme4 = fitlme(tb, 'obs4 ~ 1 + period + (1|id)');

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