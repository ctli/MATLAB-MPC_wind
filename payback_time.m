%% payback time for wind farm with battery energy storage

clear all
close all
clc
format compact


%% ================================================================== %%
data = importdata('SITE_3939_MAIN_FORECASTS.csv', ',');
obs = data.data(:,3);

scale = 800/max(obs);
obs = obs*scale;
obs = obs(1:8760)';

c1 = 1;
full_potential = -c1*sum(obs);

% energy_price = 28.96; % [$/MWh]
energy_price = 28.96 + 23; % [$/MWh]

% energy_price = 28+3+16; % [EUR/MWh], Denmark subsidy (Munksgaard  2008)


%% ================================================================== %%
overnight_wind = 1630; % [$/kW]
overnight_battery = 600; % [$/kW]

% overnight_wind = 800; % [EUR/kW]
% overnight_battery = 450; % [$/kW]

% w/o battery
total_cost0 = overnight_wind*1000*800
% annual_revenue = full_potential*0.201*energy_price; % michigan
annual_revenue = full_potential*0.36*energy_price; % Denmark
payback_period_no_battery = total_cost0/abs(annual_revenue)

disp('==============================');
% w/ Lithium battery
total_cost = overnight_wind*1000*800 + overnight_battery*1000*800*0.5
total_cost/total_cost0
% annual_revenue = full_potential*0.809*energy_price; % michigan; eta=0.95
annual_revenue = full_potential*0.8521*energy_price; % Denmark; eta=0.95
payback_period_w_battery = total_cost/abs(annual_revenue)

disp('==============================');
% w/ Sulfer battery
total_cost = overnight_wind*1000*800 + 300*1000*800*2
total_cost/total_cost0
annual_revenue = full_potential*0.7892*energy_price; % Denmark; eta=0.75
payback_period_w_battery = total_cost/abs(annual_revenue)


%% ================================================================== %%
% Denmark subsidy for new wind-power installation on land
disp('==============================');
annual_income = 2326*(2.8+0.3+1.6)/100; % [$/kW]
construction_cost = 1630; % [$/kW]
payback_denmark = construction_cost/annual_income


