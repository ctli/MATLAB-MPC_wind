%% plot x_opt & u_opt

clear all
close all
clc


%% ===================================================================== %%
data = importdata('SITE_3939_MAIN_FORECASTS.csv', ',');
obs = data.data(:,3);  %obs = reshape(obs, 24, length(obs)/24);
fcst = data.data(:,4); %fcst = reshape(fcst, 24, length(fcst)/24);

obs = obs(1:8760)';
fcst = fcst(1:8760)';

scale = 800/max(fcst);
fcst = fcst*scale;
obs = obs*scale;
err = fcst - obs; % >0:shortage, <0:excess

hr = 1:length(fcst);


%% ================================================================== %%
load ConvReserve

figure(1); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, DTE_scheduling, 'color', 'k', 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [4.2188    1.3750    3.5000    1.8542]);
set(gca, 'units', 'pixels', 'pos', [49.0144   32.5238  270.9840  134.4762]);

set(gca, 'fontsize', 7);
x_offset = 504;
% x_offset = 187;
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
ylim([-5 805]);
set(gca,'layer', 'bottom');
xlabel('Time (hr)', 'fontsize', 8);
ylabel('Wind Scheduling (MW)', 'fontsize', 8);

box off;
set(gca, 'tickdir', 'out');

[~, childObjs] = legend([h3, h1, h2], 'Wind Scheduling w/ Conventional Reserves', 'Forecast (w_f)', 'Actual (w_a)');
set(legend,'units','pixels');
set(legend, 'box', 'off');
set(legend, 'pos', [53.3333  125.0000   98.6667   53.3333]);

% export_tiff('MPC_fig7_conv_reserve', 600);


%% ================================================================== %%
cd('c_rate');
load persistence_x1_eta075
cd ..;

figure(2); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr(1:8755), u_sim, 'color', [1 0 0], 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [7.8958    1.3750    3.5000    1.8542]);
set(gca, 'units', 'pixels', 'pos', [49.0144   32.5238  270.9840  134.4762]);

set(gca, 'fontsize', 7);
x_offset = 504;
% x_offset = 187;
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
ylim([-5 805]);
set(gca,'layer', 'bottom');
xlabel('Time (hr)', 'fontsize', 8);
ylabel('Wind Scheduling (MW)', 'fontsize', 8);

box off;
set(gca, 'tickdir', 'out');

[~, childObjs] = legend([h3, h1, h2], 'Wind Scheduling w/ BESS Using the Heuristic Algorithm', 'Forecast (w_f)', 'Actual (w_a)');
set(legend,'units','pixels');
set(legend, 'box', 'off');
set(legend, 'pos', [53.3333  125.0000   120.6667   53.3333]);

% export_tiff('MPC_fig8_persistance', 600);


%% ================================================================== %%
cd('c_rate');
load Q100_N5_eta75_C4_0
cd ..;

figure(3); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr(), u_sim, 'color', [0 100 255]/255, 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    1.8542]);
set(gca, 'units', 'pixels', 'pos', [49.0144   32.5238  270.9840  134.4762]);

set(gca, 'fontsize', 7);
x_offset = 504;
% x_offset = 187;
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
ylim([-5 805]);
set(gca,'layer', 'bottom');
xlabel('Time (hr)', 'fontsize', 8);
ylabel('Wind Scheduling (MW)', 'fontsize', 8);

box off;
set(gca, 'tickdir', 'out');

[~, childObjs] = legend([h3, h1, h2], 'Wind Scheduling w/ BESS Controlled by MPC', 'Forecast (w_f)', 'Actual (w_a)');
set(legend,'units','pixels');
set(legend, 'box', 'off');
set(legend, 'pos', [53.3333  125.0000   105.6667   53.3333]);

% export_tiff('MPC_fig10_mpc', 600);

%%
figure(31); clf;
plot(hr, xn_1_sim);

set(gcf, 'unit', 'inch', 'pos', [11.5728    4.2396    3.5000    1.8542]);
set(gca, 'units', 'pixels', 'pos', [49.0144   32.5238  270.9840  134.4762]);

set(gca, 'fontsize', 7);
x_offset = 504;
% x_offset = 187;
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);

%% ==
cd('c_rate');
load Q75_N5_eta95_C4_0
cd ..;

figure(32); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr(), u_sim, 'color', [0 100 255]/255, 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    1.8542]);
set(gca, 'units', 'pixels', 'pos', [49.0144   32.5238  270.9840  134.4762]);

set(gca, 'fontsize', 7);
x_offset = 504;
% x_offset = 187;
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
ylim([-5 805]);
set(gca,'layer', 'bottom');
xlabel('Time (hr)', 'fontsize', 8);
ylabel('Wind Scheduling (MW)', 'fontsize', 8);

box off;
set(gca, 'tickdir', 'out');

[~, childObjs] = legend([h3, h1, h2], 'Wind Scheduling w/ BESS Controlled by MPC', 'Forecast (w_f)', 'Actual (w_a)');
set(legend,'units','pixels');
set(legend, 'box', 'off');
set(legend, 'pos', [53.3333  125.0000   105.6667   53.3333]);

% export_tiff('MPC_095_c4_0', 600);
figure(33); clf;
plot(hr, xn_1_sim);

set(gcf, 'unit', 'inch', 'pos', [15.2498    1.3750    3.5000    1.8542]);
set(gca, 'units', 'pixels', 'pos', [49.0144   32.5238  270.9840  134.4762]);

set(gca, 'fontsize', 7);
x_offset = 504;
% x_offset = 187;
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);


%% ================================================================== %%
cd('c_rate');
load Q75_N5_eta75_C4_55
cd ..;

figure(4); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, u_sim, 'color', [0 0.6 0], 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [15.2498    1.3750    3.5000    1.8542]);
set(gca, 'units', 'pixels', 'pos', [49.0144   32.5238  270.9840  134.4762]);

set(gca, 'fontsize', 7);
x_offset = 504;
% x_offset = 187;
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
ylim([-5 805]);
set(gca,'layer', 'bottom');
xlabel('Time (hr)', 'fontsize', 8);
ylabel('Wind Scheduling (MW)', 'fontsize', 8);

box off;
set(gca, 'tickdir', 'out');

[~, childObjs] = legend([h3, h1, h2], 'Wind Scheduling w/ BESS Controlled by Revised MPC', 'Forecast (w_f)', 'Actual (w_a)');
set(legend,'units','pixels');
set(legend, 'box', 'off');
set(legend, 'pos', [53.3333  125.0000   120.6667   53.3333]);

% export_tiff('MPC_fig12_mpc_ramp', 600);

%% ==
cd('c_rate');
load Q75_N5_eta95_C4_55
cd ..;

figure(41); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, u_sim, 'color', [0 0.6 0], 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [15.2498    1.3750    3.5000    1.8542]);
set(gca, 'units', 'pixels', 'pos', [49.0144   32.5238  270.9840  134.4762]);

set(gca, 'fontsize', 7);
x_offset = 504;
% x_offset = 187;
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
ylim([-5 805]);
set(gca,'layer', 'bottom');
xlabel('Time (hr)', 'fontsize', 8);
ylabel('Wind Scheduling (MW)', 'fontsize', 8);

box off;
set(gca, 'tickdir', 'out');

[~, childObjs] = legend([h3, h1, h2], 'Wind Scheduling w/ BESS Controlled by Revised MPC', 'Forecast (w_f)', 'Actual (w_a)');
set(legend,'units','pixels');
set(legend, 'box', 'off');
set(legend, 'pos', [53.3333  125.0000   120.6667   53.3333]);

% export_tiff('MPC_095_c4_55', 600);

