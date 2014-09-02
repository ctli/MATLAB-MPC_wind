clear all
close all
clc
format compact


%% ================================================================== %%
data = importdata('SITE_3939_MAIN_FORECASTS.csv', ',');
obs = data.data(:,3);
fcst = data.data(:,4);

scale = 800/max(fcst);
obs = obs*scale;
fcst = fcst*scale;

obs = obs(1:8760)';
fcst = fcst(1:8760)';

name_plate_capacity = 800;

std_o = std(obs);
std_o/name_plate_capacity;

dstd_o = std(diff(obs));
dstd_o/name_plate_capacity;


%% ================================================================== %%
load ConvReserve

std_DTE = std(DTE_scheduling);
std_DTE/name_plate_capacity;

dstd_DTE = std(diff(DTE_scheduling));
dstd_DTE/name_plate_capacity;


%% ================================================================== %%
file_name = {'Q75_N5_eta95_C4_0', ...
             'Q75_N5_eta95_C4_1', ...
             'Q75_N5_eta95_C4_5', ...
             'Q75_N5_eta95_C4_10', ...
             'Q75_N5_eta95_C4_20', ...
             'Q75_N5_eta95_C4_30', ...
             'Q75_N5_eta95_C4_40', ...
             'Q75_N5_eta95_C4_50', ...
             'Q75_N5_eta95_C4_55'};

curtail_pctg_table = zeros(1, length(file_name));
J1_table = zeros(1, length(file_name));
J2_table = zeros(1, length(file_name));
std_table = zeros(1, length(file_name));
dstd_table = zeros(1, length(file_name));
for i = 1:length(file_name)
    load(file_name{i});
    std_table(i) = std(u_sim);
    
    du = diff(u_sim);
    dstd_table(i) = std(du);
end

figure(1); clf;
bar(1, std_o/name_plate_capacity, 0.6, 'facec', [0 0.9 0], 'edge', 'none'); hold on;
bar(2, std_DTE/name_plate_capacity, 0.6, 'facec', [0.6 0.6 0.6], 'edge', 'none');
bar(2+(1:length(file_name)), std_table/name_plate_capacity, 0.6, 'edge', 'none');
xlim([0.5 11.5]);
ylim([0 0.302]);
set(gca, 'fontsize', 8);
set(gca, 'xtick', 1:11, 'xticklabel', {'obs', 'C.R.', '0', '0.01', '0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '0.55'});
ylabel('std of u');
xlabel('C4 range');
title('\eta=0.95');
defaultratio;
set(gcf, 'pos', [7.0729    5.8021    4.0000    3.0000]);

figure(2); clf;
bar(1, dstd_o/name_plate_capacity, 0.6, 'facec', [0 0.9 0], 'edge', 'none'); hold on;
bar(2, dstd_DTE/name_plate_capacity, 0.6, 'facec', [0.6 0.6 0.6], 'edge', 'none');
bar(2+(1:length(file_name)), dstd_table/name_plate_capacity, 0.6, 'edge', 'none');
xlim([0.5 11.5]);
ylim([0 0.302]);
set(gca, 'fontsize', 8);
set(gca, 'xtick', 1:11, 'xticklabel', {'obs', 'C.R.', '0', '0.01', '0.05', '0.1', '0.2', '0.3', '0.4', '0.5', '0.55'});
ylabel('std of du');
xlabel('C4 range');
title('\eta=0.95');
defaultratio;
set(gcf, 'pos', [11.2604    5.8021    4.0000    3.0000]);


