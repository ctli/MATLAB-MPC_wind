% control battery SOC
%
% sign convention for battery power
%   discharge: +
%   charge: -
%
% computation bottleneck: step 1 (not step 2 or 3)
%
% v3: - revise step 1, to skip repeated pre-table generation
%     - try shorter prediction horizen
% v4: - impose ramping penalty, c4*|u(k) - u(k+1)|
% v5: - different ramping penalty, c4*|u(k) - u(k-1)|
% v6: - imperfect battery efficiency, eta = 0.95 (lithium battery) & 0.75 (Sodium battery)

clear all
close all
clc
format compact


%% ================================================================== %%
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
bty_size = 0.75;%[0.25 0.5 0.75 1 2 5];
for bb = 1:length(bty_size)
battery_sizing = bty_size(bb);
battery_capacity = wind_nameplate*battery_sizing; % [MWh]
soc_range = 0:0.01:1;

% eta = 1; % battery efficiency
eta = sqrt(0.95); c_rate = 5; % battery efficiency
% eta = sqrt(0.75); c_rate = 0.15; % battery efficiency

P_dis_c_rate = battery_capacity*c_rate;
P_chg_c_rate = -P_dis_c_rate;

% sign convention for battery power
% discharge: +, charge: -
P_dis_max = battery_capacity*soc_range*eta;% when SOC is too low, discharge limit will be smaller than 800MW
P_dis_max(P_dis_max>wind_nameplate) = wind_nameplate;
P_dis_max(P_dis_max>P_dis_c_rate) = P_dis_c_rate;

P_chg_max = -battery_capacity*(1-soc_range)/eta;% when SOC is too high, charge limit will be smaller than -800MW
P_chg_max(P_chg_max<-wind_nameplate) = -wind_nameplate;
P_chg_max(P_chg_max<P_chg_c_rate) = P_chg_c_rate;


% ==========
N = 5;

du = 20;
u = 800:-du:0; u_dir = 'descend'; % (wind scheduling) better results!!
x = 0:0.02:1; % (battery SOC)

% %% coarse grid for debugging
% du = 50;
% u = 800:-du:0; u_dir = 'descend'; % (wind scheduling) better results!!
% x = 0:du/battery_capacity:1; % dx = 0.0625

cN_saturate = 128000 + 1;

c1 = 1;
c2 = 1.03;
c3 = 1;

% c4 = 0;
c4 = 0.3; % 0.5455
% c4 = 30/55; % 0.5455
c4_name = ['C4_', num2str(c4*100, '%3.0f')];

cN = 0;
% cN = cN_saturate;

x0 = 0.5;
u0 = 800;

% ==========
M = 8755; x_limit = [0 9000]; x_tick = 0:1000:9000; % 1-year horizon
% M = 7*24; x_limit = [0 24*7]; x_tick = 0:24:7*24; % 7-day horizon
% M = 72; x_limit = [0 M]; x_tick = 0:6:M; % 7-day horizon
% M = 24; x_limit = [0 24]; x_tick = 0:2:24; % 7-day horizon


% ==============================
% step 1: pre-table (instantaneous cost)
wf_range = 0:1:800;

xn_1_new = zeros(length(x), length(u), length(u), length(wf_range)); % [x1]x[x2]x[u]x[t]
J_inst_new = zeros(length(x), length(u), length(u), length(wf_range)); % [x1]x[x2]x[u]x[t]

tic;
for b = 1:length(wf_range) % loop through all wf levels
    b
    for i = 1:length(x) % loop for x1 state
        p_dis_max = interp1(soc_range, P_dis_max, x(i));
        p_dis_max = roundn(p_dis_max, -6);
        
        p_chg_max = interp1(soc_range, P_chg_max, x(i));
        p_chg_max = roundn(p_chg_max, -6);
        
        for j = 1:length(u) % loop for u control
            Rw_rqd = interp2(bin_wind_power, bin_wind_power, up_regulation_rqd, u(j), wf_range(b)); % cover shortage
            Rw_rqd = roundn(Rw_rqd, -6);
            if (Rw_rqd>p_dis_max)
                reserve_scheduling = Rw_rqd - p_dis_max;
            else
                reserve_scheduling = 0;
            end
            
            p_batt_dis = u(j) - wf_range(b); % discharge
            if (p_batt_dis>p_dis_max)
                reserve_dispatch = p_batt_dis - p_dis_max;
                curtailment = 0;
                p_batt_dis = p_dis_max;
            elseif (p_batt_dis<p_chg_max)
                reserve_dispatch = 0;
                curtailment = p_chg_max - p_batt_dis;
                p_batt_dis = p_chg_max;
            else
                reserve_dispatch = 0;
                curtailment = 0;
            end
            
            if p_batt_dis >=0
                x_new = x(i) - 1/eta*p_batt_dis/battery_capacity;
            else % p_batt_dis < 0
                x_new = x(i) - eta*p_batt_dis/battery_capacity;
            end
            x_new = roundn(x_new, -6);
            
            for ii = 1:length(u) % loop for x2 state
                J_inst_new(i,ii,j,b) = -c1*u(j) + c2*reserve_scheduling + c3*reserve_dispatch + c4*abs(u(j) - u(ii));
                xn_1_new(i,ii,j,b) = x_new;
            end % loop ii (state, x2)

        end % loop j (control, u)
    end % loop i (state, x1)
end % loop b (wf range)
toc;


% ==============================
J_sim = zeros(1,8760);
u_sim = zeros(1,8760);
xn_1_sim = zeros(1,8760);
reserve_scheduling_sim = zeros(1,8760);
reserve_dispatch1_sim = zeros(1,8760);
reserve_dispatch2_sim = zeros(1,8760);
curtailment1_sim = zeros(1,8760);
curtailment2_sim = zeros(1,8760);
p_batt = zeros(1, 8760);

x_old = x0;
u_old = u0;
tic;
for m = 1:M % control horizn
    if mod(m, 100) == 0
        m
        toc;
    end
    
    wf = fcst((1:N)+m-1);
    wa = obs((1:N)+m-1);
    
    wf = round(wf);
    
    % ==============================
    % step 2: cost table
    id_wf = interp1(wf_range, 1:length(wf_range), wf, 'nearest');
    J_inst = J_inst_new(:,:,:,id_wf);
    xn_1 = xn_1_new(:,:,:,id_wf);
    
    xn_2_temp = repmat(u, length(x), 1);
    xn_2 = zeros(size(xn_1));
    for j = 1:length(u)
        xn_2(:,:,j) = xn_2_temp;
    end
    
    J_cum = zeros(length(x), length(u), N+1); % [x1]x[x2]x[t], treat u(t-1) as another state
    u_opt = zeros(length(x), length(u), N); % [x1]x[x2]x[t]
    xn_1_opt = zeros(length(x), length(u), N); % [x1]x[x2]x[t]
    id_opt = zeros(length(x), length(u), N); % [x1]x[x2]x[t]
    
    for i = 1:length(x)
        J_cum(i,:,N+1) = cN*(x(i) - 0.5)^2;
    end
    
    for t = N:-1:1
        J_cum_temp = interp2(u,x,J_cum(:,:,t+1),xn_2(:,:,:,t),xn_1(:,:,:,t));
        
        % ensure that interplation does not give NAN
        J_cum_temp(isnan(J_cum_temp)) = cN*2;
        if any(isnan(J_cum_temp))
            disp('interpolation issue');
            pause;
        end
        
        cost = J_inst(:,:,:,t) + J_cum_temp;
        cost = roundn(cost, -6);
        [value, id] = min(cost, [], 3);
        J_cum(:,:,t) = value;
        u_opt(:,:,t) = u(id);
        id_opt(:,:,t) = id;
    end
    
    % ==============================
    % step 3: retrive trajactory
    % MPC: only impliment 1st-step control is implimented
    p_dis_max_sim = interp1(soc_range, P_dis_max, x_old);
    p_chg_max_sim = interp1(soc_range, P_chg_max, x_old);
    u_sim(m) = interp2(u,x, u_opt(:,:,1), u_old, x_old);
    J_sim(m) = interp2(u,x,J_cum(:,:,1), u_old, x_old);

    p_dis_max_sim = roundn(p_dis_max_sim, -6);
    p_chg_max_sim = roundn(p_chg_max_sim, -6);
    u_sim(m) = roundn(u_sim(m), -6);
    J_sim(m) = roundn(J_sim(m), -6);
    
    if isnan(u_sim(m))
        m
        u_sim(m)
        disp('nan issue!');
        pause;
    end
    
    Rw_rqd_sim = interp2(bin_wind_power, bin_wind_power, up_regulation_rqd, u_sim(m), wf(1)); % cover shortage
    Rw_rqd_sim = roundn(Rw_rqd_sim, -6);
    if (Rw_rqd_sim>p_dis_max_sim)
        reserve_scheduling_sim(m) = Rw_rqd_sim - p_dis_max_sim;
    else
        reserve_scheduling_sim(m) = 0;
    end
    
    % use actual wind, wa!!!
    %p_batt = u_sim(m) - wf(1); % discharge
    p_batt_dis = u_sim(m) - wa(1); % discharge
    if (p_batt_dis>p_dis_max_sim)
        reserve_dispatch1_sim(m) = p_batt_dis - p_dis_max_sim;
        p_batt_dis = p_dis_max_sim;
    elseif (p_batt_dis<p_chg_max_sim)
        curtailment1_sim(m) = p_chg_max_sim - p_batt_dis;
        p_batt_dis = p_chg_max_sim;
    end
    
    if p_batt_dis >=0
        x_new = x_old - 1/eta*p_batt_dis/battery_capacity;
    else % p_batt_dis < 0
        x_new = x_old - eta*p_batt_dis/battery_capacity;
    end
    x_new = roundn(x_new, -6);

    p_batt(t) = p_batt_dis;

    
    xn_1_sim(m) = x_new;
    x_old = x_new;
    u_old = u_sim(m);
end % loop for m (control receding)
toc;


discrepancy = sum(u_sim - p_batt + curtailment1_sim - reserve_dispatch1_sim) - sum(obs)
discrepancy_pctg = discrepancy/sum(obs)


file_name = ['Q', num2str(battery_sizing*100), '_N', num2str(N), '_eta', num2str(eta^2*100), '_', c4_name]
cd('c_rate');
save(file_name, 'm', 'battery_capacity', 'cN', 'c4', 'N', 'eta', ...
                'xn_1_sim', 'u_sim', 'J_sim', 'reserve_scheduling_sim', 'reserve_dispatch1_sim', 'reserve_dispatch2_sim', 'curtailment1_sim', 'curtailment2_sim');
cd ..;
end
            
curtail_pctg = sum(curtailment1_sim)/sum(obs)
reserve_MWh = sum(reserve_dispatch1_sim)

disp('========================================');
delta_u = diff([u0, u_sim]);

J_a = -c1*sum(u_sim)
J_1 = -c1*sum(u_sim) + c2*sum(reserve_scheduling_sim) + c3*sum(reserve_dispatch1_sim)
J_2 = -c1*sum(u_sim) + c2*sum(reserve_scheduling_sim) + c3*sum(reserve_dispatch1_sim) + sum(c4*abs(delta_u))

