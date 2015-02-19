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

battery_capacity_range = [0.25, 0.5, 0.75, 1, 2, 5, 10, 20];
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

% std_o = std(obs);
% std_o/name_plate_capacity;
% 
% dstd_o = std(diff(obs));
% dstd_o/name_plate_capacity;


%% ================================================================== %%
% load ConvReserve
load ConvReserve
total_cost_ConvReserve = sum(-c1*DTE_scheduling + c2*rw_range + c3*rw_dispatched);
temp = obs-DTE_scheduling;
temp(temp<0) = 0;
curtail_convReserve = sum(temp);
curtail_convReserve_pctg = curtail_convReserve/sum(obs);

std_DTE = std(DTE_scheduling);
std_DTE/name_plate_capacity;

dstd_DTE = std(diff(DTE_scheduling));
dstd_DTE/name_plate_capacity;


%% ================================================================== %%
J1_collection = zeros(3, 6);
J2_collection = zeros(3, 6);
curtail_pctg_collection = zeros(3, 6);
log_collection = zeros(3, 6, 3);
cd('c_rate');
for e = 1:4 % loop for eta
    switch e
        case 1
            eta = 0.75;
            file_name = {'Q25_N5_eta75_C4_55', ...
                'Q50_N5_eta75_C4_55', ...
                'Q75_N5_eta75_C4_55', ...
                'Q100_N5_eta75_C4_55', ...
                'Q200_N5_eta75_C4_55', ...
                'Q500_N5_eta75_C4_55'};
        case 2
            eta = 0.95;
            file_name = {'Q25_N5_eta95_C4_55', ...
                'Q50_N5_eta95_C4_55', ...
                'Q75_N5_eta95_C4_55', ...
                'Q100_N5_eta95_C4_55', ...
                'Q200_N5_eta95_C4_55', ...
                'Q500_N5_eta95_C4_55'};
        case 3
            eta = 0.75;
            file_name = {'Q25_N5_eta75_C4_0', ...
                'Q50_N5_eta75_C4_0', ...
                'Q75_N5_eta75_C4_0', ...
                'Q100_N5_eta75_C4_0', ...
                'Q200_N5_eta75_C4_0', ...
                'Q500_N5_eta75_C4_0'};
        case 4
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
    for i = 1:length(file_name)
        load(file_name{i});
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
    J2_collection(e,:) = J2_table;
    curtail_pctg_collection(e,:) = curtail_pctg_table;
end

%%
eta_name = {'\eta = 0.75', '\eta = 0.95'};
log_modified = log_collection;
% log_modified(3,4,2) = log_collection(3,4,3); log_modified(1,4,2) = log_modified(1,4,2) + 1.3128e+05;
% log_modified(3,4,3) = log_collection(3,4,2); log_modified(1,4,3) = log_modified(1,4,3) - 1.3128e+05;
for e = 4%1:4
    figure(e); clf; hold on;
    y = linspace(0,sum(obs)/1e6,6);
    line([0.5 7.5], [1 1]*y(2), 'color', [0.9 0.9 0.9]);
    line([0.5 7.5], [1 1]*y(3), 'color', [0.9 0.9 0.9]);
    line([0.5 7.5], [1 1]*y(4), 'color', [0.9 0.9 0.9]);
    line([0.5 7.5], [1 1]*y(5), 'color', [0.9 0.9 0.9]);
    line([0.5 7.5], [1 1]*y(6), 'color', [0.9 0.9 0.9]);
    bar(1, sum(obs)/1e6, 0.5, 'facec', [255 255 102]/255, 'edge', 'none'); hold on;
    bar(1, (1-0.62)*sum(obs)/1e6, 0.5, 'facec', [0 128 102]/255, 'edge', 'none'); hold on;
    
    b = bar(2:7, (log_modified(:,:,e))'/1e6, 0.5, 'stacked', 'edge', 'none');
    colormap summer
    
    set(gcf, 'unit', 'inch', 'pos', [6.9063+3.8*(e-1)    6.6563    3.5000    1.8500]);
    set(gca, 'pos', [0.1429    0.1771    0.8065    0.7438]);
    
    xlim([0.5 7.5]);
    set(gca, 'fontsize', 7, 'layer', 'top', 'tickdir', 'out', 'box', 'off');
    set(gca, 'xtick', 1:7, 'xticklabel', {'Conv.', '200', '400', '600', '800', '1600', '4000'});
    text(1, -0.4, 'Reserve', 'fontsize', 7, 'horizontalalignment', 'center');
    xlabel('              Battery Capacity (MWh)', 'fontsize', 8, 'pos', [3.9742   -0.33])
    yl = ylabel({'Annual Utilization of','Wind Power (Normalized)'}, 'fontsize', 8);
    
    ylim([0 sum(obs)/1e6]);
    set(gca, 'ytick', linspace(0,sum(obs)/1e6,6), 'yticklabel', linspace(0,sum(obs)/1e6,6)/(sum(obs)/1e6));

    [~, childObjs] = legend(fliplr(b), 'Curtailment', 'Battery Loss', 'Sales to Grid', 4);
    set(legend, 'units', 'pixels');
    set(legend, 'pos', [242.3333   39.2667   82.6667   48.7763]);
    lineObjs = findobj(childObjs, 'Type', 'patch');
    dx = 0.15;
    for lineIdx = 1:length(lineObjs)
        xCoords = get(lineObjs, 'XData');
        if (length(xCoords{lineIdx}) == 5)
            xlen=xCoords{lineIdx};
            set(lineObjs(lineIdx), 'XData', [xlen(1) xlen(2) xlen(3)-dx xlen(4)-dx xlen(5)]);
        end
    end
    for textIdx = 1:length(lineObjs)
        textObjs = findobj(childObjs, 'Type', 'text');
        textpos=get(textObjs(textIdx),'position');
        set(textObjs(textIdx), 'position', [textpos(1)-dx, textpos(2:3)]);
    end

        if e==1
        text(1, (whole-(0.62)*sum(obs)/2)/1e6, '62%', 'fontsize', 7, 'horizontalalignment', 'center');
        text(2, (whole-log_modified(3,1,e)/2)/1e6, num2str(log_modified(3,1,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(3, (whole-log_modified(3,2,e)/2)/1e6, num2str(log_modified(3,2,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(4, (whole-log_modified(3,3,e)/2)/1e6, num2str(log_modified(3,3,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(5, (whole-log_modified(3,4,e)/2)/1e6, num2str(log_modified(3,4,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(6, (whole-log_modified(3,5,e)/2)/1e6, num2str(log_modified(3,5,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');

        text(2, (log_modified(1,1,e)+log_modified(2,1,e)/2)/1e6, num2str(log_modified(2,1,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(3, (log_modified(1,2,e)+log_modified(2,2,e)/2)/1e6, num2str(log_modified(2,2,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(4, (log_modified(1,3,e)+log_modified(2,3,e)/2)/1e6, num2str(log_modified(2,3,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(5, (log_modified(1,4,e)+log_modified(2,4,e)/2)/1e6, num2str(log_modified(2,4,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(6, (log_modified(1,5,e)+log_modified(2,5,e)/2)/1e6, num2str(log_modified(2,5,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(7, (log_modified(1,6,e)+log_modified(2,6,e)/2)/1e6, num2str(log_modified(2,6,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        
        text(7.35, 2.31, '(BESS eff: \eta=0.87)', 'fontsize', 8, 'horizontalalignment', 'right');
        % export_fig MPC_fig13_eta075_v2 -painters
        end

        if e==4
        text(1, (whole-(0.62)*sum(obs)/2)/1e6, '62%', 'fontsize', 7, 'horizontalalignment', 'center');
        text(2, (whole-log_modified(3,1,e)/2)/1e6, num2str(log_modified(3,1,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(3, (whole-log_modified(3,2,e)/2)/1e6, num2str(log_modified(3,2,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(4, (whole-log_modified(3,3,e)/2)/1e6, num2str(log_modified(3,3,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        text(5, (whole-log_modified(3,4,e)/2)/1e6, num2str(log_modified(3,4,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
%         text(6, (whole-log_modified(3,5,e)/2)/1e6, '<1%', 'fontsize', 7, 'horizontalalignment', 'center');
%         text(7, (whole-log_modified(3,6,e)/2)/1e6, '<1%', 'fontsize', 7, 'horizontalalignment', 'center');

        text(2, (log_modified(1,1,e)+log_modified(2,1,e)/2)/1e6, '<1%', 'fontsize', 7, 'horizontalalignment', 'center');
        text(3, (log_modified(1,2,e)+log_modified(2,2,e)/2)/1e6, '<1%', 'fontsize', 7, 'horizontalalignment', 'center');
        text(4, (log_modified(1,3,e)+log_modified(2,3,e)/2)/1e6, '<1%', 'fontsize', 7, 'horizontalalignment', 'center');
        text(5, (log_modified(1,4,e)+log_modified(2,4,e)/2)/1e6, num2str(log_modified(2,4,e)/whole*100, '%2.0f%%'), 'fontsize', 7, 'horizontalalignment', 'center');
        
        text(7.35, 2.31, '(BESS eff: \eta=0.97)', 'fontsize', 8, 'horizontalalignment', 'right');
        % export_fig MPC_fig13_eta095_v2 -painters
        end
end

