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

wa = obs;
wf = fcst;

wind_nameplate = 800; % [MW]

hr = 1:length(fcst);

%% ================================================================== %%
battery_capacity_range = [0.25, 0.5, 0.75, 1, 2, 5, 10];
soc_range = 0:0.05:1;

J1_collection = zeros(2, length(battery_capacity_range));
J2_collection = zeros(2, length(battery_capacity_range));

tic;
for e = 1:2 % loop for {eta=0.75, eta=0.95}

    if e==1
        eta = sqrt(0.75); c_rate = 0.15; eta_name = 'eta075'; alpha = 0.8;
    end
    
    if e==2
        eta = sqrt(0.95); c_rate = 5; eta_name = 'eta095'; alpha = 0.9;
    end
    
    J1_table = zeros(1, 6);
    J2_table = zeros(1, 6);
    
for q = 1:length(battery_capacity_range)

    battery_capacity = wind_nameplate*battery_capacity_range(q)
    P_dis_c_rate = battery_capacity*c_rate;
    P_chg_c_rate = -P_dis_c_rate;

    P_dis_max = battery_capacity*soc_range; % when SOC is too low, discharge limit will be smaller than 800MW
    P_dis_max(P_dis_max>wind_nameplate) = wind_nameplate;
    P_dis_max(P_dis_max>P_dis_c_rate) = P_dis_c_rate;

    P_chg_max = -battery_capacity*(1-soc_range); % when SOC is too high, charge limit will be smaller than -800MW
    P_chg_max(P_chg_max<-wind_nameplate) = -wind_nameplate;
    P_chg_max(P_chg_max<P_chg_c_rate) = P_chg_c_rate;
    
    soc_ini = 0.5;
    xn_1_sim = zeros(1,8755);
    reserve_scheduling_sim = zeros(1,8755);
    reserve_dispatch1_sim = zeros(1,8755);
    curtailment1_sim = zeros(1,8755);
    
    x0 = soc_ini;
    x_old = x0;
    
    tmp1 = [wa(1), wa(1:8755-1)];
    tmp2 = wf(1:8755);
    u_sim = alpha*tmp1 + (1-alpha)*tmp2;
    Rw_rqd_sim = interp2(bin_wind_power, bin_wind_power, up_regulation_rqd, u_sim, wf(1:8755)); % cover shortage
    
    for m = 1:8755 % loop for time
        p_dis_max_sim = interp1(soc_range, P_dis_max, x_old);
        p_dis_max_sim = roundn(p_dis_max_sim, -6);
        
        p_chg_max_sim = interp1(soc_range, P_chg_max, x_old);
        p_chg_max_sim = roundn(p_chg_max_sim, -6);
        
        if (Rw_rqd_sim(m)>p_dis_max_sim)
            reserve_scheduling_sim(m) = Rw_rqd_sim(m) - p_dis_max_sim;
        else
            reserve_scheduling_sim(m) = 0;
        end
        
        p_batt_dis = u_sim(m) - wa(m); % discharge
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
        x_new = min(x_new,1);
        x_new = max(x_new,0);
        
        xn_1_sim(m) = x_new;
        x_old = x_new;
    end % loof of m
    
    curtail_pctg = sum(curtailment1_sim)/sum(obs);
    reserve_MWh = sum(reserve_dispatch1_sim);
    
    c1 = 1;
    c2 = 1.03;
    c3 = 1;
    cost = -c1*u_sim + c2*reserve_scheduling_sim + c3*reserve_dispatch1_sim;
    J1 = sum(cost)

    c4 = 30/55;
    u0 = 400;
    delta_u = diff([u0, u_sim]);
    cost_corrected = -c1*u_sim + c2*reserve_scheduling_sim + c3*reserve_dispatch1_sim + c4*abs(delta_u);
    J2 = sum(cost_corrected);

    J1_table(q) = J1;
    J2_table(q) = J2;
    
    cd('c_rate');
    if battery_capacity_range(q) < 1
        temp_name = ['0', num2str(battery_capacity_range(q)*100)];
    else
        temp_name = num2str(battery_capacity_range(q));
    end
    file_name = ['persistence_x', temp_name, '_', eta_name];
    save(file_name, 'battery_capacity_range', 'xn_1_sim', 'u_sim', 'reserve_scheduling_sim', 'reserve_dispatch1_sim', 'curtailment1_sim');
    cd ..;
    toc;
    disp('==============================');
end % loop of q, battery capacity

    J1_collection(e,:) = J1_table;
    J2_collection(e,:) = J2_table;
    
end % loop of e, battery efficiency

full_potential = -c1*sum(obs);

%% ================================================================== %%
load ConvReserve
total_cost_ConvReserve = sum(-c1*DTE_scheduling + c2*rw_range + c3*rw_dispatched);
temp = obs-DTE_scheduling;
temp(temp<0) = 0;
curtail_convReserve = sum(temp);
curtail_convReserve_pctg = curtail_convReserve/sum(obs);
% ==============================

figure(1); clf; hold on;
h1 = bar(1, total_cost_ConvReserve/full_potential, 0.5, 'facec', [0.3 0.3 0.3], 'edge', 'none'); hold on;
h2 = bar(2:length(J1_collection)+1, J1_collection'/full_potential, 1, 'group', 'edge', 'none');
set(h2(1), 'facec', [1 0 0]);
set(h2(2), 'facec', [255 110 110]/255);

set(gcf, 'unit', 'inch', 'pos', [6.9063    0.8646    3.5000    1.8500]);
set(gca, 'pos', [0.1429    0.1771    0.8065    0.7555]);
xlim([0.5 7.5]);
ylim([0 1.01]);
set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
ylabel('Annual Revenue (Normalized)', 'fontsize', 8);
set(gca, 'xtick', 1:7, 'xticklabel', {'Conv.', '200', '400', '600', '800', '1600', '4000'});
text(1, -0.16, 'Reserve', 'fontsize', 7, 'horizontalalignment', 'center');
xlabel('              Battery Capacity (MWh)', 'fontsize', 8, 'pos', [3.9742   -0.15]);

