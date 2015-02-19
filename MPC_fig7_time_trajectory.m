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
load Q25_N5_eta75_C4_0
cd ..;

figure(3); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, u_sim, 'color', [0 100 255]/255, 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    1.8542]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.1761    0.8065    0.7513]);

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

% export_tiff('MPC_fig10_mpc_a', 600);

% ==========
figure(31); clf;
plot(hr, xn_1_sim, '-', 'color', [0 100 255]/255, 'linewidth', 1);

box off;
set(gca, 'fontsize', 7);
ylim([-0.02 1.02]);
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
xlabel('Time (hr)', 'fontsize', 8);
ylabel('SOC (-)', 'fontsize', 8);

set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    1.8542]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.1761    0.8065    0.7513]);

% ==========
figure(32); clf;
plot(xn_1_sim, u_sim/800, '.');

set(gca, 'fontsize', 7);
axis equal;
axis([0 1 0 1]);
set(gca, 'xtick', 0:0.2:1, 'ytick', 0:0.2:1);
xlabel('SOC (-)', 'fontsize', 8);
ylabel('u (-)', 'fontsize', 8);
title(corr2(xn_1_sim, u_sim/800));
set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    2.25]);

%% ==============================
cd('c_rate');
load Q500_N5_eta75_C4_0
u1 = u_sim;

load Q200_N5_eta75_C4_0
u2 = u_sim;
cd ..;

aa = 0.9;
u3 = aa*u1 + (1-aa)*u2;

figure(4); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, u3, 'color', [0 100 255]/255, 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [11.5728    4.2188    3.5000    1.8542]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.1761    0.8065    0.7513]);

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

% export_tiff('MPC_fig10_mpc_b', 600);

% ==========
cd('c_rate');
load Q500_N5_eta75_C4_0
cd ..;

figure(41); clf;
plot(hr, xn_1_sim, '-', 'color', [0 100 255]/255, 'linewidth', 1);

box off;
set(gca, 'fontsize', 7);
ylim([-0.02 1.02]);
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
xlabel('Time (hr)', 'fontsize', 8);
ylabel('SOC (-)', 'fontsize', 8);

set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    1.8542]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.1761    0.8065    0.7513]);

% ==========

figure(42); clf;
plot(xn_1_sim, u_sim/800, '.');

set(gca, 'fontsize', 7);
axis equal;
axis([0 1 0 1]);
set(gca, 'xtick', 0:0.2:1, 'ytick', 0:0.2:1);
xlabel('SOC (-)', 'fontsize', 8);
ylabel('u (-)', 'fontsize', 8);
title(corr2(xn_1_sim, u_sim/800));
set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    2.25]);

%% ==============================
cd('c_rate');
load Q100_N5_eta75_C4_0
u1 = u_sim;

load Q200_N5_eta75_C4_0
u2 = u_sim;
cd ..;

u3 = 0.5*u1+0.5*u2;

figure(5); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, u3, 'color', [0 100 255]/255, 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [11.5728    4.2188    3.5000    1.8542]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.1761    0.8065    0.7513]);

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

% export_tiff('MPC_fig10_mpc_bb', 600);

% ==========
cd('c_rate');
load Q200_N5_eta75_C4_0
cd ..;

figure(51); clf;
plot(hr, xn_1_sim, '-', 'color', [0 100 255]/255, 'linewidth', 1);

box off;
set(gca, 'fontsize', 7);
ylim([-0.02 1.02]);
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
xlabel('Time (hr)', 'fontsize', 8);
ylabel('SOC (-)', 'fontsize', 8);

set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    1.8542]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.1761    0.8065    0.7513]);

% ==========
figure(52); clf;
plot(xn_1_sim, u_sim/800, '.');

set(gca, 'fontsize', 7);
axis equal;
axis([0 1 0 1]);
set(gca, 'xtick', 0:0.2:1, 'ytick', 0:0.2:1);
xlabel('SOC (-)', 'fontsize', 8);
ylabel('u (-)', 'fontsize', 8);
title(corr2(xn_1_sim, u_sim/800));
set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    2.25]);

%% ==============================
cd('c_rate');
load Q25_N5_eta95_C4_0
cd ..;

figure(6); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, u_sim, 'color', [0 100 255]/255, 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [11.5728    7.0625    3.5000    1.8542]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.1761    0.8065    0.7513]);

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

% export_tiff('MPC_fig10_mpc_c', 600);

% ==========
figure(61); clf;
plot(hr, xn_1_sim, '-', 'color', [0 100 255]/255, 'linewidth', 1);

box off;
set(gca, 'fontsize', 7);
ylim([-0.02 1.02]);
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
xlabel('Time (hr)', 'fontsize', 8);
ylabel('SOC (-)', 'fontsize', 8);

set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    1.8542]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.1761    0.8065    0.7513]);

% ==========
figure(62); clf;
plot(xn_1_sim, u_sim/800, '.');

set(gca, 'fontsize', 7);
axis equal;
axis([0 1 0 1]);
set(gca, 'xtick', 0:0.2:1, 'ytick', 0:0.2:1);
xlabel('SOC (-)', 'fontsize', 8);
ylabel('u (-)', 'fontsize', 8);
title(corr2(xn_1_sim, u_sim/800));
set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    2.25]);

%% ==============================
cd('c_rate');
load Q200_N5_eta95_C4_0
cd ..;

figure(7); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, u_sim, 'color', [0 100 255]/255, 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [11.5728    7.0625    3.5000    1.8542]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.1761    0.8065    0.7513]);

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

% export_tiff('MPC_fig10_mpc_d', 600);

% ==========
figure(71); clf;
plot(hr, xn_1_sim, '-', 'color', [0 100 255]/255, 'linewidth', 1);

box off;
set(gca, 'fontsize', 7);
ylim([-0.02 1.02]);
xlim([0 7*24+0.3]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
xlabel('Time (hr)', 'fontsize', 8);
ylabel('SOC (-)', 'fontsize', 8);

set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    1.8542]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.1761    0.8065    0.7513]);

% ==========
figure(72); clf;
plot(xn_1_sim, u_sim/800, '.');

set(gca, 'fontsize', 7);
axis equal;
axis([0 1 0 1]);
set(gca, 'xtick', 0:0.2:1, 'ytick', 0:0.2:1);
xlabel('SOC (-)', 'fontsize', 8);
ylabel('u (-)', 'fontsize', 8);
title(corr2(xn_1_sim, u_sim/800));
set(gcf, 'unit', 'inch', 'pos', [11.5728    1.3750    3.5000    2.25]);

%% ================================================================== %%
cd('c_rate');
load Q25_N5_eta75_C4_55
cd ..;

figure(8); clf;
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

% export_tiff('MPC_fig12_mpc_ramp_a', 600);

%% ====
cd('c_rate');
load Q50_N5_eta95_C4_55
cd ..;

figure(9); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, u_sim, 'color', [0 0.6 0], 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [15.2498    4.2292    3.5000    1.8542]);
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

% export_tiff('MPC_fig12_mpc_ramp_b', 600);


%% ====
cd('c_rate');
load Q75_N5_eta95_C4_55
cd ..;

figure(10); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, u_sim, 'color', [0 0.6 0], 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [15.2498    4.2292    3.5000    1.8542]);
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

% export_tiff('MPC_fig12_mpc_ramp_b', 600);

% ====
cd('c_rate');
load Q75_N5_eta95_C4_30
cd ..;

figure(11); clf;
h1 = plot(hr, fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
h2 = plot(hr, obs, '--', 'color', [0.1 0.1 0.1]);
h3 = plot(hr, u_sim, 'color', [0 0.6 0], 'linewidth', 1);

set(gcf, 'unit', 'inch', 'pos', [15.2498    4.2292    3.5000    1.8542]);
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

% export_tiff('MPC_fig12_mpc_ramp_30', 600);
