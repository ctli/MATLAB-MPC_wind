clear all
close all
clc

wind_nameplate = 800; % [MW]

bty_size = [0.25 0.5 0.75 1 2 5];

bb = 5%1:length(bty_size)
battery_sizing = bty_size(bb);
battery_capacity = wind_nameplate*battery_sizing; % [MWh]
soc_range = 0:0.01:1;

% eta = 1; % battery efficiency
% eta = sqrt(0.95); % battery efficiency
eta = sqrt(0.75); % battery efficiency

% sign convention for battery power
% discharge: +, charge: -
P_dis_max = battery_capacity*soc_range*eta;% when SOC is too low, discharge limit will be smaller than 800MW
P_dis_max(P_dis_max>wind_nameplate) = wind_nameplate;

P_chg_max = -battery_capacity*(1-soc_range)/eta;% when SOC is too high, charge limit will be smaller than -800MW
P_chg_max(P_chg_max<-wind_nameplate) = -wind_nameplate;


figure(1); clf;
plot(soc_range, P_dis_max);
hold on;
plot(soc_range, P_chg_max);

c_rate = 5;
P_chg_c_rate = battery_capacity*c_rate;
P_dis_c_rate = -P_chg_c_rate;
plot([0 1], [1 1]*P_chg_c_rate, 'r');
plot([0 1], [1 1]*P_dis_c_rate), 'r');
