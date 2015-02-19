clear all
close all
clc
format compact


%% ================================================================== %%
c1 = 1;
c2 = 1.03;
c3 = 1;
c4 = 30/55;

u0 = 400;

name_plate_capacity = 800;


%% ================================================================== %%
data = importdata('SITE_3939_MAIN_FORECASTS.csv', ',');
obs = data.data(:,3);
fcst = data.data(:,4);

scale = 800/max(fcst);
obs = obs*scale;
fcst = fcst*scale;

obs = obs(1:8760)';
fcst = fcst(1:8760)';

wind_nameplate = 800; % [MW]

full_potential = -c1*sum(obs);


% fcst = fcst(1:8755);


%% ================================================================== %%
cd('c_rate');


%% Step change
figure(1); clf; hold on; box on;
b = -400:50:400;

% ========
load('persistence_x50_eta095');
du0 = diff(u_sim);
n0 = hist(du0, b);
figure(1); bar(b,n0/length(du0), 'facec', 'r', 'linestyle', 'none');

figure(3);
plot(xn_1_sim, u_sim, '.');
corr2(xn_1_sim/800, u_sim)

figure(31);
plot(fcst(1:length(u_sim)), u_sim, '.');
corr2(fcst(1:length(u_sim)), u_sim)

disp('====================');
load('Q200_N5_eta95_C4_55');
du1 = diff(u_sim);
n1 = hist(du1, b);
figure(1); bar(b,n1/length(du1), 0.5, 'facec', 'g', 'linestyle', 'none');

corr2(xn_1_sim/800, u_sim)
corr2(fcst(1:length(u_sim)), u_sim)

