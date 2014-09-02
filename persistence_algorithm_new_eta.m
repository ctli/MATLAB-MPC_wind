%% persistance algorithm for wind + battery

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
% dn-regulatin requirement, surplus, need charge battery, <0
dn_regulation_rqd = zeros(length(bin_wind_power), length(bin_wind_power)); % [fcst]x[scheduling] = [wf]x[u2]
for i = 1:length(bin_wind_power)
    for j = 1:length(bin_wind_power)
        u2 = bin_wind_power(j);
        temp = u2 - excess_generation(i);
        temp = min(0, temp);
        dn_regulation_rqd(i,j) = temp;
    end
end

% ==============================
% expected deviation
Rw_expected_dispatch = zeros(length(bin_wind_power), length(bin_wind_power)); % [fcst]x[scheduling] = [wf]x[u2]
for j = 1:length(bin_wind_power) % [scheduling]
    u2 = bin_wind_power(j);
    for i = 1:length(bin_wind_power) % [fcst]
        temp = u2 - bin_wind_power;
        %temp = max(temp, 0);
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

% wa = obs;
% wf = fcst;

wind_nameplate = 800; % [MW]


%% ================================================================== %%
battery_capacity_range = [0.25, 0.5, 0.75, 1, 2, 5];
soc_range = 0:0.05:1;

% eta = 0.95; eta_name = 'eta095';
eta = 0.75; eta_name = 'eta075';

tic;
for q = 1:length(battery_capacity_range)
    battery_capacity = wind_nameplate*battery_capacity_range(q)
    
    P_dis_max = battery_capacity*soc_range; % when SOC is too low, discharge limit will be smaller than 800MW
    P_dis_max(P_dis_max>wind_nameplate) = wind_nameplate;
    
    P_chg_max = -battery_capacity*(1-soc_range); % when SOC is too high, charge limit will be smaller than -800MW
    P_chg_max(P_chg_max<-wind_nameplate) = -wind_nameplate;
    
    soc_ini = 0.5;
    xn_1_sim = zeros(1,8760);
    u_sim = zeros(1,8760);
    reserve_scheduling_sim = zeros(1,8760);
    reserve_dispatch1_sim = zeros(1,8760);
    reserve_dispatch2_sim = zeros(1,8760);
    curtailment1_sim = zeros(1,8760);
    curtailment2_sim = zeros(1,8760);
    
    x0 = soc_ini;
    x_old = x0;
    
    u0 = obs(1);
    u_old = u0;
    for m = 1:8755 % loop for time
        wf = fcst(m);
        wa = obs(m);
        
        p_dis_max_sim = interp1(soc_range, P_dis_max, x_old);
        p_dis_max_sim = roundn(p_dis_max_sim, -6);
        
        p_chg_max_sim = interp1(soc_range, P_chg_max, x_old);
        p_chg_max_sim = roundn(p_chg_max_sim, -6);
        
        u = (u_old + fcst(m))/2; % persistence strategy
        
        Rw_rqd_sim = interp2(bin_wind_power, bin_wind_power, up_regulation_rqd, u, wf); % cover shortage
        Rw_rqd_sim = roundn(Rw_rqd_sim, -6);
        if (Rw_rqd_sim>p_dis_max_sim)
            reserve_scheduling_sim(m) = Rw_rqd_sim - p_dis_max_sim;
        else
            reserve_scheduling_sim(m) = 0;
        end
        
        p_batt_dis = u - wa; % discharge
        if (p_batt_dis>p_dis_max_sim)
            reserve_dispatch1_sim(m) = p_batt_dis - p_dis_max_sim;
            p_batt_dis = p_dis_max_sim;
        elseif (p_batt_dis<p_chg_max_sim)
            curtailment1_sim(m) = p_chg_max_sim - p_batt_dis;
            p_batt_dis = p_chg_max_sim;
        end
        
        if p_batt_dis >=0
            x_new = x_old - 1/eta*p_batt_dis/battery_capacity;
        else % charging
            x_new = x_old - eta*p_batt_dis/battery_capacity;
        end
        x_new = roundn(x_new, -6);
        
        u_sim(m) = u;
        xn_1_sim(m) = x_new;
        x_old = x_new;
        u_old = wa;
        
    end % loof of m
    
    curtail_pctg = sum(curtailment1_sim)/sum(obs)
    reserve_MWh = sum(reserve_dispatch1_sim)
    
    c1 = 1;
    c2 = 1.03;
    c3 = 1;
    cost = -c1*u_sim + c2*reserve_scheduling_sim + c3*reserve_dispatch1_sim;
    cost_sum = sum(cost)

    c4 = 30/55;
    u0 = 400;
    delta_u = diff([u0, u_sim]);
    cost_corrected = -c1*u_sim + c2*reserve_scheduling_sim + c3*reserve_dispatch1_sim + c4*abs(delta_u);
    cost_corrected_sum = sum(cost_corrected)
    
    if battery_capacity_range(q) < 1
        temp_name = ['0', num2str(battery_capacity_range(q)*100)];
    else
        temp_name = num2str(battery_capacity_range(q));
    end
    file_name = ['persistence_x', temp_name, '_', eta_name];
    save(file_name, 'battery_capacity_range', 'xn_1_sim', 'u_sim', 'reserve_scheduling_sim', 'reserve_dispatch1_sim', 'reserve_dispatch2_sim', 'curtailment1_sim', 'curtailment2_sim');
    
    toc;
    disp('==============================');
end % loop of q, battery capacity



%% ================================================================== %%
figure(1); clf;
plot(fcst, '-', 'color', [0.75 0.75 0.75]); hold on;
plot(obs, '--', 'color', [0.1 0.1 0.1]);
plot(u_sim, 'b');

x_offset = 504;
xlim([0 7*24]+x_offset);
set(gca, 'fontsize', 8);
xlabel('Time (hr)');
ylabel('Wind Scheduling (MW)');

