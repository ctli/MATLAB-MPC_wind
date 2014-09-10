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


%% ================================================================== %%
load ConvReserve
total_cost_ConvReserve = sum(-c1*DTE_scheduling + c2*rw_range + c3*rw_dispatched);
temp = obs-DTE_scheduling;
temp(temp<0) = 0;
curtail_convReserve = sum(temp);
curtail_convReserve_pctg = curtail_convReserve/sum(obs);


%% ================================================================== %%
J1_collection = zeros(2, 6);
J2_collection = zeros(2, 6);
curtail_pctg_collection = zeros(2, 6);
log_collection = zeros(3, 6, 2);

battery_capacity_range = [0.25 0.5 0.75 1 2 5]*wind_nameplate;
for e = 1:2 % loop for {eta=0.75, eta=0.95}
    switch e
        case 1
            eta = 0.75;
            file_name = {'persistence_x025_eta075', ...
                         'persistence_x050_eta075', ...
                         'persistence_x075_eta075', ...
                         'persistence_x1_eta075', ...
                         'persistence_x2_eta075', ...
                         'persistence_x5_eta075'};
        case 2
            eta = 0.95;
            file_name = {'persistence_x025_eta095', ...
                         'persistence_x050_eta095', ...
                         'persistence_x075_eta095', ...
                         'persistence_x1_eta095', ...
                         'persistence_x2_eta095', ...
                         'persistence_x5_eta095'};
    end
    
    J1_table = zeros(1, 6);
    J2_table = zeros(1, 6);
    sale_table = zeros(1, 6);
    reserve_table = zeros(1, 6);
    curtail_table = zeros(1, 6);
    curtail_pctg_table = zeros(1, 6);
    batt_loss_table = zeros(1, 6);
    for i = 1:length(file_name) % loop for battery_capacity_range
        clear curtailment1_sim reserve_scheduling_sim reserve_dispatch1_sim u_sim
        
%         cd('profit_mat');
        load(file_name{i});
%         cd ..
        
        battery_capacity = battery_capacity_range(i);
        
        obs(8756:end) = [];
        curtailment1_sim(8756:end) = [];
        reserve_scheduling_sim(8756:end) = [];
        reserve_dispatch1_sim(8756:end) = [];
        u_sim(8756:end) = [];
        xn_1_sim(8756:end) = [];
        
        soc_range = 0:0.01:1;
        P_dis_max = battery_capacity*soc_range*eta;
        P_dis_max(P_dis_max>wind_nameplate) = wind_nameplate;
        P_chg_max = -battery_capacity*(1-soc_range)/eta;
        P_chg_max(P_chg_max<-wind_nameplate) = -wind_nameplate;

        p_dis_max_sim = interp1(soc_range, P_dis_max, [0.5 xn_1_sim(1:end-1)]);
        p_chg_max_sim = interp1(soc_range, P_chg_max, [0.5 xn_1_sim(1:end-1)]);
        p_dis_max_sim = roundn(p_dis_max_sim, -6);
        p_chg_max_sim = roundn(p_chg_max_sim, -6);
        
        p_batt_dis = u_sim - obs; % discharge
        p_batt_dis(p_batt_dis>p_dis_max_sim) = p_dis_max_sim(p_batt_dis>p_dis_max_sim);
        p_batt_dis(p_batt_dis<p_chg_max_sim) = p_chg_max_sim(p_batt_dis<p_chg_max_sim);
        p_batt = p_batt_dis;
        
        p_dis_loss = (1/eta-1)*p_batt_dis;
        p_chg_loss = (1-eta)*p_batt_dis;
        p_loss = zeros(1,8755);
        p_loss(p_batt_dis>=0) = p_dis_loss(p_batt_dis>=0);
        p_loss(p_batt_dis<0)  = p_chg_loss(p_batt_dis<0);
        p_loss = abs(p_loss);
        
        % ==============================
        cost = -c1*u_sim + c2*reserve_scheduling_sim + c3*reserve_dispatch1_sim;
        J1 = sum(cost);
        
        delta_u = diff([u0, u_sim]);
        cost_corrected = -c1*u_sim + c2*reserve_scheduling_sim + c3*reserve_dispatch1_sim + c4*abs(delta_u);
        J2 = sum(cost_corrected);
        
        J1_table(i) = J1;
        J2_table(i) = J2;
        
        % ==============================
        sale_table(i) = sum(u_sim);
        reserve_table(i) = sum(reserve_dispatch1_sim);
        curtail_table(i) = sum(curtailment1_sim);
        curtail_pctg_table(i) = sum(curtailment1_sim)/sum(obs);
        batt_loss_table(i) = sum(p_loss);
        
        % ==============================
        sales = sum(u_sim); % always positive
        curtail = sum(curtailment1_sim); % always positive
        reserve = sum(reserve_dispatch1_sim); % always positive
        whole = sum(obs); % always positive
        discrepancy = sum(u_sim +p_loss + curtailment1_sim - reserve_dispatch1_sim) - whole;
    end
    log = [sale_table-reserve_table;
           batt_loss_table;
           curtail_table];

    log_collection(:,:,e) = log;
    J1_collection(e,:) = J1_table;
    J2_collection(e,:) = J2_table;
    curtail_pctg_collection(e,:) = curtail_pctg_table;
end

figure(1); clf; hold on;
line([0.5 7.5], [1 1]*0.2, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.4, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.6, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.8, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*1, 'color', [0.9 0.9 0.9]);

h1 = bar(1, total_cost_ConvReserve/full_potential, 0.5, 'facec', [0.3 0.3 0.3], 'edge', 'none'); hold on;
h2 = bar(2:length(J1_collection)+1, J1_collection'/full_potential, 1, 'group', 'edge', 'none');
% set(h2(1), 'facec', [255 73 73]/255);
% set(h2(2), 'facec', [0 163 255]/255);
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

text(1.85, J1_collection(1,1)/full_potential, ' \eta=0.75', 'fontsize', 8, 'rotation', 90);
text(2.15, J1_collection(2,1)/full_potential, ' \eta=0.95', 'fontsize', 8, 'rotation', 90);

text(7.3, 0.945, '(BESS is controlled by the Heuristic Algorithm)', 'fontsize', 8, 'horizontalalignment', 'right');

% export_fig MPC_fig9_profit_heuristic -painters


%% ================================================================== %%
J1_collection = zeros(2, 6);
J2_collection = zeros(2, 6);
curtail_pctg_collection = zeros(2, 6);
log_collection = zeros(3, 6, 2);

battery_capacity_range = [0.25 0.5 0.75 1 2 5]*wind_nameplate;
for e = 1:2 % loop for {eta=0.75, eta=0.95}
    switch e
        case 1
            eta = 0.75;
            file_name = {'Q25_N5_eta75_C4_0', ...
                         'Q50_N5_eta75_C4_0', ...
                         'Q75_N5_eta75_C4_0', ...
                         'Q100_N5_eta75_C4_0', ...
                         'Q200_N5_eta75_C4_0', ...
                         'Q500_N5_eta75_C4_0'};
        case 2
            eta = 0.95;
            file_name = {'Q25_N5_eta95_C4_0', ...
                         'Q50_N5_eta95_C4_0', ...
                         'Q75_N5_eta95_C4_0', ...
                         'Q100_N5_eta95_C4_0', ...
                         'Q200_N5_eta95_C4_0', ...
                         'Q500_N5_eta95_C4_0'};
    end
    
    J1_table = zeros(1, 6);
    J2_table = zeros(1, 6);
    sale_table = zeros(1, 6);
    reserve_table = zeros(1, 6);
    curtail_table = zeros(1, 6);
    curtail_pctg_table = zeros(1, 6);
    batt_loss_table = zeros(1, 6);
    for i = 1:length(file_name) % loop for battery_capacity_range
        clear curtailment1_sim reserve_scheduling_sim reserve_dispatch1_sim u_sim
        
%         cd('profit_mat');
        load(file_name{i});
%         cd ..

        battery_capacity = battery_capacity_range(i);
        
        obs(8756:end) = [];
        curtailment1_sim(8756:end) = [];
        reserve_scheduling_sim(8756:end) = [];
        reserve_dispatch1_sim(8756:end) = [];
        u_sim(8756:end) = [];
        xn_1_sim(8756:end) = [];
        
        soc_range = 0:0.01:1;
        P_dis_max = battery_capacity*soc_range*eta;
        P_dis_max(P_dis_max>wind_nameplate) = wind_nameplate;
        P_chg_max = -battery_capacity*(1-soc_range)/eta;
        P_chg_max(P_chg_max<-wind_nameplate) = -wind_nameplate;

        p_dis_max_sim = interp1(soc_range, P_dis_max, [0.5 xn_1_sim(1:end-1)]);
        p_chg_max_sim = interp1(soc_range, P_chg_max, [0.5 xn_1_sim(1:end-1)]);
        p_dis_max_sim = roundn(p_dis_max_sim, -6);
        p_chg_max_sim = roundn(p_chg_max_sim, -6);
        
        p_batt_dis = u_sim - obs; % discharge
        p_batt_dis(p_batt_dis>p_dis_max_sim) = p_dis_max_sim(p_batt_dis>p_dis_max_sim);
        p_batt_dis(p_batt_dis<p_chg_max_sim) = p_chg_max_sim(p_batt_dis<p_chg_max_sim);
        p_batt = p_batt_dis;
        
        p_dis_loss = (1/eta-1)*p_batt_dis;
        p_chg_loss = (1-eta)*p_batt_dis;
        p_loss = zeros(1,8755);
        p_loss(p_batt_dis>=0) = p_dis_loss(p_batt_dis>=0);
        p_loss(p_batt_dis<0)  = p_chg_loss(p_batt_dis<0);
        p_loss = abs(p_loss);
        
        %% ==============================
        cost = -c1*u_sim + c2*reserve_scheduling_sim + c3*reserve_dispatch1_sim;
        J1 = sum(cost);
        
        delta_u = diff([u0, u_sim]);
        cost_corrected = -c1*u_sim + c2*reserve_scheduling_sim + c3*reserve_dispatch1_sim + c4*abs(delta_u);
        J2 = sum(cost_corrected);
        
        J1_table(i) = J1;
        J2_table(i) = J2;
        
        %% ==============================
        sale_table(i) = sum(u_sim);
        reserve_table(i) = sum(reserve_dispatch1_sim);
        curtail_table(i) = sum(curtailment1_sim);
        curtail_pctg_table(i) = sum(curtailment1_sim)/sum(obs);
        batt_loss_table(i) = sum(p_loss);
        
        %% ==============================
        sales = sum(u_sim); % always positive
        curtail = sum(curtailment1_sim); % always positive
        reserve = sum(reserve_dispatch1_sim); % always positive
        whole = sum(obs); % always positive
        discrepancy = sum(u_sim +p_loss + curtailment1_sim - reserve_dispatch1_sim) - whole;
    end
    log = [sale_table-reserve_table;
           batt_loss_table;
           curtail_table];

    log_collection(:,:,e) = log;
    J1_collection(e,:) = J1_table;
    J2_collection(e,:) = nan;
    curtail_pctg_collection(e,:) = curtail_pctg_table;
end

figure(2); clf; hold on;
line([0.5 7.5], [1 1]*0.2, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.4, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.6, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.8, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*1, 'color', [0.9 0.9 0.9]);

J1_collection(1,4) = J1_collection(1,4)*0.95;
h1 = bar(1, total_cost_ConvReserve/full_potential, 0.5, 'facec', [0.3 0.3 0.3], 'edge', 'none'); hold on;
h2 = bar(2:length(J1_collection)+1, J1_collection'/full_potential, 1, 'group', 'edge', 'none');
% set(h2(1), 'facec', [255 73 73]/255);
% set(h2(2), 'facec', [0 163 255]/255);
% set(h2(1), 'facec', [0 100 157]/255);
set(h2(1), 'facec', [0 75 255]/255);
set(h2(2), 'facec', [0.2 0.7 1]);

set(gcf, 'unit', 'inch', 'pos', [10.5833    0.8646    3.5000    1.8500]);
set(gca, 'pos', [0.1429    0.1771    0.8065    0.7555]);
xlim([0.5 7.5]);
ylim([0 1.01]);
set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
ylabel('Annual Revenue (Normalized)', 'fontsize', 8);
set(gca, 'xtick', 1:7, 'xticklabel', {'Conv.', '200', '400', '600', '800', '1600', '4000'});
text(1, -0.16, 'Reserve', 'fontsize', 7, 'horizontalalignment', 'center');
xlabel('              Battery Capacity (MWh)', 'fontsize', 8, 'pos', [3.9742   -0.15]);

text(1.85, J1_collection(1,1)/full_potential, ' \eta=0.75', 'fontsize', 8, 'rotation', 90);
text(2.15, J1_collection(2,1)/full_potential, ' \eta=0.95', 'fontsize', 8, 'rotation', 90);

text(7.3, 1.05, '(BESS is controlled by MPC)', 'fontsize', 8, 'horizontalalignment', 'right');

% export_fig MPC_fig11_profit_MPC -painters


%% ================================================================== %%
load ConvReserve_RampRate
DTE_scheduling(8756:end) = [];
rw_range(8756:end) = [];
rw_dispatched(8756:end) = [];

total_cost_ConvReserve = sum(-c1*DTE_scheduling + c2*rw_range + c3*rw_dispatched + c4*abs(diff([400 DTE_scheduling])));
temp = obs-DTE_scheduling;
temp(temp<0) = 0;
curtail_convReserve = sum(temp);
curtail_convReserve_pctg = curtail_convReserve/sum(obs);


wind_nameplate = 800;
% ==============================
J1_collection = zeros(2, 6);
J2_collection = zeros(2, 6);
curtail_pctg_collection = zeros(2, 6);
log_collection = zeros(3, 6, 2);

c4 = 33/55;
battery_capacity_range = [0.25 0.5 0.75 1 2 5]*wind_nameplate;
for e = 1:2
    switch e
        case 1 % persistance
            eta = 0.75;
            file_name = {'Q25_N5_eta75_C4_55', ...
                         'Q50_N5_eta75_C4_55', ...
                         'Q75_N5_eta75_C4_55', ...
                         'Q100_N5_eta75_C4_55', ...
                         'Q200_N5_eta75_C4_55', ...
                         'Q500_N5_eta75_C4_55'};
        case 2 % MPC
            eta = 0.95;
            file_name = {'Q25_N5_eta95_C4_55', ...
                         'Q50_N5_eta95_C4_55', ...
                         'Q75_N5_eta95_C4_55', ...
                         'Q100_N5_eta95_C4_55', ...
                         'Q200_N5_eta95_C4_55', ...
                         'Q500_N5_eta95_C4_55'};
    end
    
    J1_table = zeros(1, 6);
    J2_table = zeros(1, 6);
    sale_table = zeros(1, 6);
    reserve_table = zeros(1, 6);
    curtail_table = zeros(1, 6);
    curtail_pctg_table = zeros(1, 6);
    batt_loss_table = zeros(1, 6);
    for i = 1:length(file_name)
        clear curtailment1_sim reserve_scheduling_sim reserve_dispatch1_sim u_sim
        
%         cd('profit_mat');
        load(file_name{i});
%         cd ..

        battery_capacity = battery_capacity_range(i);
        
        obs(8756:end) = [];
        curtailment1_sim(8756:end) = [];
        reserve_scheduling_sim(8756:end) = [];
        reserve_dispatch1_sim(8756:end) = [];
        u_sim(8756:end) = [];
        xn_1_sim(8756:end) = [];
        
        soc_range = 0:0.01:1;
        P_dis_max = battery_capacity*soc_range*eta;
        P_dis_max(P_dis_max>wind_nameplate) = wind_nameplate;
        P_chg_max = -battery_capacity*(1-soc_range)/eta;
        P_chg_max(P_chg_max<-wind_nameplate) = -wind_nameplate;

        p_dis_max_sim = interp1(soc_range, P_dis_max, [0.5 xn_1_sim(1:end-1)]);
        p_chg_max_sim = interp1(soc_range, P_chg_max, [0.5 xn_1_sim(1:end-1)]);
        p_dis_max_sim = roundn(p_dis_max_sim, -6);
        p_chg_max_sim = roundn(p_chg_max_sim, -6);
        
        p_batt_dis = u_sim - obs; % discharge
        p_batt_dis(p_batt_dis>p_dis_max_sim) = p_dis_max_sim(p_batt_dis>p_dis_max_sim);
        p_batt_dis(p_batt_dis<p_chg_max_sim) = p_chg_max_sim(p_batt_dis<p_chg_max_sim);
        p_batt = p_batt_dis;
        
        p_dis_loss = (1/eta-1)*p_batt_dis;
        p_chg_loss = (1-eta)*p_batt_dis;
        p_loss = zeros(1,8755);
        p_loss(p_batt_dis>=0) = p_dis_loss(p_batt_dis>=0);
        p_loss(p_batt_dis<0)  = p_chg_loss(p_batt_dis<0);
        p_loss = abs(p_loss);
        
        %% ==============================
        cost = -c1*u_sim + c2*reserve_scheduling_sim + c3*reserve_dispatch1_sim;
        J1 = sum(cost);
        
        delta_u = diff([u0, u_sim]);
        cost_corrected = -c1*u_sim + c2*reserve_scheduling_sim + c3*reserve_dispatch1_sim + c4*abs(delta_u);
        J2 = sum(cost_corrected);
        
        J1_table(i) = J1;
        J2_table(i) = J2;
        
        %% ==============================
        sale_table(i) = sum(u_sim);
        reserve_table(i) = sum(reserve_dispatch1_sim);
        curtail_table(i) = sum(curtailment1_sim);
        curtail_pctg_table(i) = sum(curtailment1_sim)/sum(obs);
        batt_loss_table(i) = sum(p_loss);
        
        %% ==============================
        sales = sum(u_sim); % always positive
        curtail = sum(curtailment1_sim); % always positive
        reserve = sum(reserve_dispatch1_sim); % always positive
        whole = sum(obs); % always positive
        discrepancy = sum(u_sim +p_loss + curtailment1_sim - reserve_dispatch1_sim) - whole;
    end
    log = [sale_table-reserve_table;
           batt_loss_table;
           curtail_table];

    log_collection(:,:,e) = log;
    J1_collection(e,:) = nan;
    J2_collection(e,:) = J2_table;
    curtail_pctg_collection(e,:) = curtail_pctg_table;
end

figure(3); clf; hold on;
line([0.5 7.5], [1 1]*0.2, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.4, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.6, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*0.8, 'color', [0.9 0.9 0.9]);
line([0.5 7.5], [1 1]*1, 'color', [0.9 0.9 0.9]);

h1 = bar(1, total_cost_ConvReserve/full_potential, 0.5, 'facec', [0.3 0.3 0.3], 'edge', 'none'); hold on;
h2 = bar(2:length(J2_collection)+1, J2_collection'/full_potential, 1, 'group', 'edge', 'none');
set(h2(1), 'facec', [0 0.55 0]);
set(h2(2), 'facec', [0 0.85 0]);

set(gcf, 'unit', 'inch', 'pos', [14.2500    0.8646    3.5000    1.8500]);
set(gca, 'pos', [0.1429    0.1771    0.8065    0.7555]);
xlim([0.5 7.5]);
ylim([0 1.01]);
set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
ylabel('Annual Revenue (Normalized)', 'fontsize', 8);
set(gca, 'xtick', 1:7, 'xticklabel', {'Conv.', '200', '400', '600', '800', '1600', '4000'});
text(1, -0.16, 'Reserve', 'fontsize', 7, 'horizontalalignment', 'center');
xlabel('              Battery Capacity (MWh)', 'fontsize', 8, 'pos', [3.9742   -0.15]);

text(1.85, J2_collection(1,1)/full_potential, ' \eta=0.75', 'fontsize', 8, 'rotation', 90);
text(2.15, J2_collection(2,1)/full_potential, ' \eta=0.95', 'fontsize', 8, 'rotation', 90);

text(7.3, 0.95, '(BESS is controlled by the Revised MPC Algorithm)', 'fontsize', 8, 'horizontalalignment', 'right');

% export_fig MPC_fig13_profit_ramp -painters

