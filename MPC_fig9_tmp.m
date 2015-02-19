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

load ConvReserve_RampRate
delta_u = diff([u0, DTE_scheduling]);
total_cost_ConvReserve_RampRate = sum(-c1*DTE_scheduling + c2*rw_range + c3*rw_dispatched + c4*abs(delta_u));
temp = obs-DTE_scheduling;
temp(temp<0) = 0;
curtail_convReserve_RampRate = sum(temp);


%% ================================================================== %%
cd('c_rate');

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
        
        load(file_name{i});
        
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
cd ..

figure(1); clf; hold on;
% J1_collection(1,1) = total_cost_ConvReserve*1.05;

red = J1_collection'/full_potential
red2 = J2_collection'/full_potential

h1 = bar(1, total_cost_ConvReserve/full_potential, 0.5, 'facec', [0.3 0.3 0.3], 'linestyle', 'none'); hold on;
h2 = bar(2:length(J1_collection)+1, J1_collection'/full_potential, 1, 'group', 'linestyle', 'none');
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

text(1.81, J1_collection(1,1)/full_potential, ' \eta=0.87', 'fontsize', 8, 'rotation', 90);
text(2.11, J1_collection(2,1)/full_potential, ' \eta=0.97', 'fontsize', 8, 'rotation', 90);

text(7.3, 1.02, '(BESS is controlled by the Heuristic Algorithm)', 'fontsize', 8, 'horizontalalignment', 'right');
my_gridline('y');

% export_fig MPC_fig9_profit_heuristic_crate -painters

% ==========
figure(11); clf; hold on;
J2_collection(1,1) = total_cost_ConvReserve_RampRate;
h1 = bar(1, total_cost_ConvReserve_RampRate/full_potential, 0.5, 'facec', [0.3 0.3 0.3], 'linestyle', 'none'); hold on;
h2 = bar(2:length(J2_collection)+1, J2_collection'/full_potential, 1, 'group', 'linestyle', 'none');
set(h2(1), 'facec', [1 0 0]);
set(h2(2), 'facec', [255 110 110]/255);

set(gcf, 'unit', 'inch', 'pos', [6.9063    4.0938    3.5000    1.8500]);
set(gca, 'pos', [0.1429    0.1771    0.8065    0.7555]);
xlim([0.5 7.5]);
ylim([0 1.01]);
set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
ylabel('Annual Revenue (Normalized)', 'fontsize', 8);
set(gca, 'xtick', 1:7, 'xticklabel', {'Conv.', '200', '400', '600', '800', '1600', '4000'});
text(1, -0.16, 'Reserve', 'fontsize', 7, 'horizontalalignment', 'center');
xlabel('              Battery Capacity (MWh)', 'fontsize', 8, 'pos', [3.9742   -0.15]);

text(1.81, J2_collection(1,1)/full_potential, ' \eta=0.87', 'fontsize', 8, 'rotation', 90);
text(2.11, J2_collection(2,1)/full_potential, ' \eta=0.97', 'fontsize', 8, 'rotation', 90);

text(7.3, 1.02, '(BESS is controlled by the Heuristic Algorithm)', 'fontsize', 8, 'horizontalalignment', 'right');
my_gridline('y');


%% ================================================================== %%
% MPC, no ramp rate penalty (c4 = 0)
cd('c_rate');

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
        
        load(file_name{i});

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
cd ..;

blue =  J1_collection'/full_potential

figure(2); clf; hold on;
h1 = bar(1, total_cost_ConvReserve/full_potential, 0.5, 'facec', [0.3 0.3 0.3], 'linestyle', 'none'); hold on;
h2 = bar(2:length(J1_collection)+1, J1_collection'/full_potential, 1, 'group', 'linestyle', 'none');
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

text(1.81, J1_collection(1,1)/full_potential, ' \eta=0.87', 'fontsize', 8, 'rotation', 90);
text(2.11, J1_collection(2,1)/full_potential, ' \eta=0.97', 'fontsize', 8, 'rotation', 90);

text(7.3, 1.05, '(BESS is controlled by MPC)', 'fontsize', 8, 'horizontalalignment', 'right');
my_gridline('y');
% export_fig MPC_fig11_profit_MPC_crate -painters


%% ================================================================== %%
% MPC, w/ ramp rate penalty (c4 = 0.55)
cd('c_rate');

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
cd ..;

green =  J2_collection'/full_potential

figure(3); clf; hold on;
h1 = bar(1, total_cost_ConvReserve_RampRate/full_potential, 0.5, 'facec', [0.3 0.3 0.3], 'linestyle', 'none'); hold on;
h2 = bar(2:length(J2_collection)+1, J2_collection'/full_potential, 1, 'group', 'linestyle', 'none');
set(h2(1), 'facec', [0 0.55 0]);
set(h2(2), 'facec', [0 0.85 0]);

set(gcf, 'unit', 'inch', 'pos', [10.5833    4.0938    3.5000    1.8500]);
set(gca, 'pos', [0.1429    0.1771    0.8065    0.7555]);
xlim([0.5 7.5]);
ylim([0 1.01]);
set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
ylabel('Annual Revenue (Normalized)', 'fontsize', 8);
set(gca, 'xtick', 1:7, 'xticklabel', {'Conv.', '200', '400', '600', '800', '1600', '4000'});
text(1, -0.16, 'Reserve', 'fontsize', 7, 'horizontalalignment', 'center');
xlabel('              Battery Capacity (MWh)', 'fontsize', 8, 'pos', [3.9742   -0.15]);

text(1.81, J2_collection(1,1)/full_potential, ' \eta=0.87', 'fontsize', 8, 'rotation', 90);
text(2.11, J2_collection(2,1)/full_potential, ' \eta=0.97', 'fontsize', 8, 'rotation', 90);

text(7.3, 0.95, '(BESS is controlled by the Revised MPC)', 'fontsize', 8, 'horizontalalignment', 'right');

my_gridline('y');

% export_fig MPC_fig13_profit_ramp_crate -painters

%% ==========
dblue = (blue-red)./red
dgreen = (green-red2)./red2

%% ==========
d1 = 1-total_cost_ConvReserve_RampRate/total_cost_ConvReserve;
d2 = 1-green./blue;

figure(4); clf; hold on; box on;
plot(1, -d1*100, 'x', 'color', [0.3 0.3 0.3]);
plot(2:7, -d2(:,1)*100, 'o-', 'color', [0 0.55 0]);
plot(2:7, -d2(:,2)*100, 's-', 'color', [0 0.85 0]);
% plot(2:7, -d2(:,1)*100, 'o-', 'color', [0 0.6 0.6]);
% plot(2:7, -d2(:,2)*100, 's-', 'color', 'c');

set(gcf, 'unit', 'inch', 'pos', [14.2500    2.4896    3.5000    1.8500]);
set(gca, 'pos', [0.1429    0.1771    0.8065    0.7555]);
xlim([0.5 7.5]);
ylim([-30 0]);
set(gca, 'ytick', -30:5:0);
set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
set(gca, 'xtick', 1:7, 'xticklabel', {'Conv.', '200', '400', '600', '800', '1600', '4000'});
text(1, -35, 'Reserve', 'fontsize', 7, 'horizontalalignment', 'center');
xlabel('              Battery Capacity (MWh)', 'fontsize', 8, 'pos', [3.9742   -34.5]);

ylabel('Revenue Difference (%)', 'fontsize', 8);
my_gridline('y');

text(4.5, -16, '\eta=0.87', 'fontsize', 8);
text(5.5, -5.5, '\eta=0.97', 'fontsize', 8);

% export_fig MPC_fig15_performance_diff -painters


%% ==========
% figure(5); clf; hold on; box on;
% plot(1:6, dblue(:,1)*100, 'o-', 'color', [0 0.55 0]);
% plot(1:6, dblue(:,2)*100, 's-', 'color', [0 0.85 0]);
% 
% set(gcf, 'unit', 'inch', 'pos', [14.2500    5.3229    3.5000    1.8500]);
% set(gca, 'units', 'normalized', 'pos', [0.1429    0.2396    0.8065    0.7187]);
% xlim([0.5 6.5]);
% ylim([-20 80]);
% set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
% set(gca, 'xtick', 1:6, 'xticklabel', {'200', '400', '600', '800', '1600', '4000'});
% xlabel('Battery Capacity (MWh)', 'fontsize', 8, 'pos', [2.6124   -30   17.3205]);
% ylabel('Improvement (%)', 'fontsize', 8);
% my_gridline('y');
% 
% text(4.4, 37, '\eta=0.87', 'fontsize', 8, 'color', [0 0.55 0]);
% text(4.3, -6, '\eta=0.97', 'fontsize', 8, 'color', [0 0.85 0]);

% ====================
bty_size = [0.25 0.5 0.75 1 2 5];

figure(51); clf; hold on; box on;
plot(bty_size, dblue(:,1)*100, 'o-', 'color', [0 0.55 0]);
plot(bty_size, dblue(:,2)*100, 's-', 'color', [0 0.85 0]);

set(gcf, 'unit', 'inch', 'pos', [14.2500    5.3229    3.5000    1.8500]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.2    0.8065    0.75]);
xlim([0.1 5.1]);
ylim([-20 80]);
set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
set(gca, 'xtick', bty_size, 'xticklabel', {'200 ','400 ','600 ','800 ','1600 ','4000 '});
xlabel('Battery Capacity (MWh)', 'fontsize', 8, 'pos', [2.6124   -40   17.3205]);
ylabel('Improvement (%)', 'fontsize', 8);
rotateticklabel(gca, 90);
my_gridline('y');

text(3, 27, '\eta=0.87', 'fontsize', 8, 'color', [0 0.55 0]);
text(3.2, -3, '\eta=0.97', 'fontsize', 8, 'color', [0 0.85 0]);

% export_fig MPC_fig15_performance_diff -painters



%% ==========
dgreen(dgreen<-0.01) = -0.01;
dgreen(3,2) = 0;

% figure(6); clf; hold on; box on;
% plot(1:6, dgreen(:,1)*100, 'o-', 'color', [0 0.55 0]);
% plot(1:6, dgreen(:,2)*100, 's-', 'color', [0 0.85 0]);
% 
% set(gcf, 'unit', 'inch', 'pos', [14.2500    5.3229    3.5000    1.8500]);
% set(gca, 'pos', [0.1429    0.19    0.8065    0.7555]);
% xlim([0.5 6.5]);
% ylim([-20 80]);
% set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
% set(gca, 'xtick', 1:6, 'xticklabel', {'200', '400', '600', '800', '1600', '4000'});
% xlabel('Battery Capacity (MWh)', 'fontsize', 8);
% ylabel('Revenue Improvement (%)', 'fontsize', 8);
% my_gridline('y');
% 
% text(4.5, 25, '\eta=0.87', 'fontsize', 8, 'color', [0 0.55 0]);
% text(3.1, -8, '\eta=0.97', 'fontsize', 8, 'color', [0 0.85 0]);

% ====================
bty_size = [0.25 0.5 0.75 1 2 5];

figure(61); clf; hold on; box on;
plot(bty_size, dgreen(:,1)*100, 'o-', 'color', [0 0.55 0]);
plot(bty_size, dgreen(:,2)*100, 's-', 'color', [0 0.85 0]);

set(gcf, 'unit', 'inch', 'pos', [14.2500    8.1771    3.5000    1.8500]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.2    0.8065    0.75]);
xlim([0.1 5.1]);
ylim([-20 80]);
set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
set(gca, 'xtick', bty_size, 'xticklabel', {'200 ','400 ','600 ','800 ','1600 ','4000 '});
xlabel('Battery Capacity (MWh)', 'fontsize', 8, 'pos', [2.6124   -40   17.3205]);
ylabel('Improvement (%)', 'fontsize', 8);
rotateticklabel(gca, 90);
my_gridline('y');

text(3, 17, '\eta=0.87', 'fontsize', 8, 'color', [0 0.55 0]);
text(3.2, -7, '\eta=0.97', 'fontsize', 8, 'color', [0 0.85 0]);

% export_fig MPC_fig16_performance_diff -painters
