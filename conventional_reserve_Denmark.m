clear all
close all
clc
format compact


%% ================================================================== %%
% wind energy & reserve requirment
bin_size = 50;
bin_wind_power = 0:bin_size:800;

load('NREL_probability_map'); % 'group_cdf', 'group_pdf', 'group_mean', 'prob_trans_matrix'
load('guarantee_generation'); % 'guarantee_generation'

% ==============================
% up-regulatin requirement, shortage, need discharge battery, >0
up_regulation_rqd = zeros(length(bin_wind_power), length(bin_wind_power)); % [fcst]x[scheduling] = [wf]x[u2]
for i = 1:length(bin_wind_power)
    for j = 1:length(bin_wind_power)
        u2 = bin_wind_power(j);
        temp = u2 - guarantee_generation(i);
        temp = max(0, temp);
        up_regulation_rqd(i,j) = temp;
    end
end

% ==============================
% expected reserve dispatch
Rw_expected_dispatch = zeros(length(bin_wind_power), length(bin_wind_power)); % [fcst]x[scheduling] = [wf]x[u2]
for j = 1:length(bin_wind_power) % [scheduling]
    u2 = bin_wind_power(j);
    for i = 1:length(bin_wind_power) % [fcst]
        temp = u2 - bin_wind_power;
        temp = max(temp, 0); % no battery
        temp_e = sum(temp.*group_pdf(i,:));
        Rw_expected_dispatch(i,j) = temp_e;
    end
end

% ==============================
data = importdata('SITE_3939_MAIN_FORECASTS.csv', ',');
obs = data.data(:,3);
fcst = data.data(:,4);

scale = 800/max(fcst);
obs = obs*scale;
fcst = fcst*scale;

obs = obs(1:8760)';
fcst = fcst(1:8760)';


%% ================================================================== %%
% wind scheduling optimization
c1 = 28+3+16;
c2 = 28*1.03;
c3 = 28*1;

% Denmark wind subsidy:
% J. Munksgaard and P. Morthorst, "Wind power in the Danish liberalised 
% power market - policy measures, price impact and investor incentives,
% " Energy Policy, vol. 36, pp. 3940-3947, 2008.

scheduling_decision = zeros(1, length(bin_wind_power)); % [1]x[fcst]
opt_cost_map = zeros(1, length(bin_wind_power)); % [1]x[fcst]
Rw_rqd_map = zeros(1, length(bin_wind_power)); % [1]x[fcst]
rw_dispatched_map = zeros(1, length(bin_wind_power)); % [1]x[fcst]
for f = 1:length(bin_wind_power) % loop for [fcst]
    r_schedule = up_regulation_rqd(f,:); % 1x[scheduling]
    r_disaptched = Rw_expected_dispatch(f,:); % 1x[scheduling]
    
    total_cost = -c1*bin_wind_power + c2*r_schedule + c3*r_disaptched;
    [value, id] = min(total_cost);
    opt_cost_map(f) = value;
    scheduling_decision(f) = bin_wind_power(id);
    Rw_rqd_map(f) = r_schedule(id);
    rw_dispatched_map(f) = r_disaptched(id);
end % end of loop f (loop for [fcst])

DTE_scheduling = interp1(bin_wind_power, scheduling_decision, fcst);
rw_range = interp1(bin_wind_power, Rw_rqd_map, fcst);

rw_dispatched = DTE_scheduling-obs;
rw_dispatched(rw_dispatched<0) = 0;

total_cost = sum(-c1*DTE_scheduling + c2*rw_range + c3*rw_dispatched)

% save('ConvReserve_Denmark', 'DTE_scheduling', 'rw_range', 'rw_dispatched');

% ==============================
hr = 1:length(fcst);
figure(1); clf;
plot(hr, fcst, hr, obs, hr, DTE_scheduling);

set(gcf, 'unit', 'inch', 'pos', [8.3229    5.1771    3.4000    2.0000]);
set(gca, 'units', 'pixels', 'pos', [47.9899   33.1824  271.0176  140.1038]);
set(gca, 'fontsize', 8);
x_offset = 504;
% x_offset = 0;
xlim([0 7*24]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
ylim([-20 820]);
xlabel('Time (hr)');
ylabel('u*, Wind Scheduling (MW)');


%% ================================================================== %%
% wind scheduling optimization, with ramp rate penalty
scheduling_decision = zeros(length(bin_wind_power), length(bin_wind_power)); % [uk]x[fcst]
opt_cost_map = zeros(length(bin_wind_power), length(bin_wind_power)); % [uk]x[fcst]
Pl_reduced_map = zeros(length(bin_wind_power), length(bin_wind_power)); % [uk]x[fcst]
Rw_rqd_map = zeros(length(bin_wind_power), length(bin_wind_power)); % [uk]x[fcst]
rw_dispatched_map = zeros(length(bin_wind_power), length(bin_wind_power)); % [uk]x[fcst]

c4 = 28*0.25;
for uk = 1:length(bin_wind_power) % loop for [uk]
    u_previous = bin_wind_power(uk);
    for f = 1:length(bin_wind_power) % loop for [fcst]
        r_schedule = up_regulation_rqd(f,:); % 1x[scheduling]
        r_disaptched = Rw_expected_dispatch(f,:); % 1x[scheduling]
        
        total_cost = -c1*bin_wind_power + c2*r_schedule + c3*r_disaptched + c4*abs(bin_wind_power-u_previous);
        [value id] = min(total_cost);
        opt_cost_map(uk,f) = value;
        scheduling_decision(uk,f) = bin_wind_power(id);
        Rw_rqd_map(uk,f) = r_schedule(id);
        rw_dispatched_map(uk,f) = r_disaptched(id);
    end % end of loop f (loop for [fcst])
end % end of loop uk (loop for [uk])

u0 = fcst(1);
DTE_scheduling = zeros(1, length(fcst));
rw_range = zeros(1, length(fcst));
uk = u0;
for t = 1:length(fcst)
    DTE_scheduling(t) = interp2(bin_wind_power, bin_wind_power, scheduling_decision, fcst(t), uk);
    rw_range(t) = interp2(bin_wind_power, bin_wind_power, Rw_rqd_map, fcst(t), uk);
    uk = DTE_scheduling(t);
end
rw_dispatched = DTE_scheduling-obs;
rw_dispatched(rw_dispatched<0) = 0;

delta_u = diff([u0, DTE_scheduling]);
total_cost = sum(-c1*DTE_scheduling + c2*rw_range + c3*rw_dispatched + c4*abs(delta_u))

full_potential = -c1*sum(obs)

ratio = total_cost/full_potential

% save('ConvReserve_Denmark_RampRate', 'DTE_scheduling', 'rw_range', 'rw_dispatched', 'c4', 'u0');

% ==============================
hr = 1:length(fcst);
figure(2); clf;
plot(hr, fcst, hr, obs, hr, DTE_scheduling);

set(gcf, 'unit', 'inch', 'pos', [8.3229    2.1979    3.4000    2.0000]);
set(gca, 'units', 'pixels', 'pos', [47.9899   33.1824  271.0176  140.1038]);
set(gca, 'fontsize', 8);
x_offset = 504;
% x_offset = 0;
xlim([0 7*24]+x_offset);
set(gca, 'xtick', (0:24:8760));
set(gca, 'xticklabel', (0:24:8760)-x_offset);
ylim([-20 820]);
xlabel('Time (hr)');
ylabel('u*, Wind Scheduling (MW)');


