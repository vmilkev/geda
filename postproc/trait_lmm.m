
clear;

traits1 = {"trait_v1_05-Mar-2025_19-58-33_robot_101.geda",...
          "trait_v1_05-Mar-2025_19-59-09_robot_102.geda",...
          "trait_v1_05-Mar-2025_19-59-44_robot_103.geda"};
robots1 = [101, 102, 103];

traits4 = {"trait_v1_05-Mar-2025_15-31-42_robot_1.geda",...
           "trait_v1_05-Mar-2025_15-36-21_robot_2.geda",...
           "trait_v1_05-Mar-2025_15-41-47_robot_3.geda"};
robots4 = [1, 2, 3];

traits6 = {"trait_v1_05-Mar-2025_19-49-23_robot_101.geda"};
robots6 = [101];

[tr] = fit_model(traits1, robots1, 1);

disp("-----------------------------------------");

%% Nested Functions

function [tr] = fit_model(trait_list, robot_list, fig_id)

tmp_tables = cell( size(trait_list,2), 1 );
for i = 1:size(trait_list,2)
    tmp_tables{i,1} = import_trait(trait_list{i});
    robot = zeros( height(tmp_tables{i,1}), 1 );
    robot(:,1) = robot_list(i);
    tmp_tables{i,1}.("robot") = robot;
end

tr = vertcat(tmp_tables{:});

u = unique(tr.id);
obs_per_id = size(tr,1)/numel(u);

tn = convertTo(tr.time, 'posixtime');
tn = round((tn - min(tn) + 1)./60);
tc = zeros( size(tn) );
h = 8;
period = h*60; % observation period in 24 hours (n_perionds = (24/h)*60), min
day = 24*60;
for i = 1:size(tn,1)
    rel_time = floor((tn(i,1) - floor(tn(i,1)/day)*day)/period);
    tc(i,1) = rel_time + 1;
end

tday = tr.time;
tday.Format = 'dd-MMM-yyyy';
tday = categorical(string(tday));

tr.("period") = tc;
tr.("day") = tday;

% f = find( abs(tr.bkg_1(:)) > 0.3 );
% tr(f,:) = [];

if size(robot_list,2) > 1
    model1 = 'gas_1 ~ 1 + period + day + robot * bkg_1 * period + (1|id)';    
    model2 = 'gas_2 ~ 1 + period + day + robot * bkg_2 * period + (1|id)';
else
    model1 = 'gas_1 ~ 1 + period + day + bkg_1 + (1|id)';    
    model2 = 'gas_2 ~ 1 + period + day + bkg_2 + (1|id)';
end

lme = fitlme(tr, model1);
[psi,mse] = covarianceParameters(lme);
hco4 = psi{1}/(psi{1}+mse); % repeatability/heritability
disp(["sigmma", "err", "h", "r^2", "n_obs_id"]);
disp([psi{1}, mse, psi{1}/(psi{1}+mse), lme.Rsquared.Adjusted, obs_per_id]);

lme2 = fitlme(tr, model2);
[psi,mse] = covarianceParameters(lme2);
%disp([psi{1}, mse, psi{1}/(psi{1}+mse), lme2.Rsquared.Adjusted]);

if fig_id == 0
    return;
end

figure(fig_id);
clf;
set(gcf,'Color','white');
% --- prediction CH4 ----------
subplot(2,2,1)
iobs = tr.gas_1;
yhat = predict(lme,tr);
x = min(iobs):max(iobs);
y = min(iobs):max(iobs);
hold on;
plot(iobs, yhat, '*');
plot(x,y,'-','LineWidth',1.5);
xlim([min(iobs) max(iobs)]);
ylim([min(iobs) max(iobs)]);
xlabel("CH4 observation");
ylabel("CH4 prediction");
box on
axis square
title( strcat("CH4: adjusted r^2 = ", num2str(lme.Rsquared.Adjusted)) );
% --- prediction CO2 ----------
subplot(2,2,2)
iobs = tr.gas_2;
yhat = predict(lme2,tr);
x = min(iobs):max(iobs);
y = min(iobs):max(iobs);
hold on;
plot(iobs, yhat, '*');
plot(x,y,'-','LineWidth',1.5);
xlim([min(iobs) max(iobs)]);
ylim([min(iobs) max(iobs)]);
xlabel("CO2 observation");
ylabel("CO2 prediction");
box on
axis square
title( strcat("CO2: adjusted r^2 = ", num2str(lme2.Rsquared.Adjusted)) );
% --- correlation CH4 vs bkg ----------
subplot(2,2,3)
x = linspace(min(tr.gas_1),max(tr.gas_1),10);
y = linspace(max(tr.bkg_1),min(tr.bkg_1),10);
hold on
plot(tr.gas_1,tr.bkg_1, '*');
% plot(x, y, 'r-','LineWidth',1.5);
% xlim([min(x) max(x)]);
% ylim([min(y) max(y)]);
xlabel("CH4 trait");
ylabel("CH4 bkg corection val.");
box on
axis square
title( strcat("CH4: corr = ", num2str(corr(tr.gas_1,tr.bkg_1))) );
% --- correlation CO2 vs bkg ----------
subplot(2,2,4)
plot(tr.gas_2,tr.bkg_2, '*')
xlabel("CO2 trait");
ylabel("CO2 bkg corection val.");
box on
axis square
title( strcat("CO2: corr = ", num2str(corr(tr.gas_2,tr.bkg_2))) );
sgtitle( strcat( "LMM: estimated repeatability h^2_{CH4} = ", num2str(hco4) ) );

end

function tr1 = import_trait(filename, dataLines)
if nargin < 2
    dataLines = [2, Inf];
end
opts = delimitedTextImportOptions("NumVariables", 6);
opts.DataLines = dataLines;
opts.Delimiter = ";";
opts.VariableNames = ["time", "id", "gas_1", "bkg_1", "gas_2", "bkg_2"];
opts.VariableTypes = ["datetime", "uint64", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "time", "InputFormat", "dd-MMM-yyyy HH:mm:ss");
tr1 = readtable(filename, opts);
end