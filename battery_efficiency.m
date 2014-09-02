%% test battery efficiency

clear all
close all
clc
format compact

Battery_int; % Prius battery

figure(1); clf;
plot(ess_soc, ess_max_pwr_dis/1e3, ...
     ess_soc, ess_max_pwr_chg/1e3);
xlabel('SOC (-)');
ylabel('Power (kW)');
legend('Discharge', 'Charge', 4);
title('Discharge/Charge Power Limits');
defaultratio;
set(gcf, 'pos', [2.4583    5.7708    4.0000    3.0000]);


%% ================================================================== %%
Pdis = (1:30)*1e3; % [kW], terminal power, i.e. output power
[SOC,P_DIS] = meshgrid(ess_soc, Pdis);
VOC = interp1(ess_soc, ess_voc, SOC);
R_DIS = interp1(ess_soc, ess_r_dis, SOC);

I_DIS = (VOC - sqrt(VOC.^2 - 4*P_DIS.*R_DIS))./(2*R_DIS);
I_DIS(abs(imag(I_DIS))>0) = nan;
Ptotal = VOC.*I_DIS;
% Ploss = (I_DIS.^2).*R_DIS;

eff_dis = P_DIS./Ptotal;
% eff_dis = P_DIS./(P_DIS+Ploss);

ave_eff_dis = nanmean(eff_dis(:))

figure(2); clf; 
mesh(ess_soc, Pdis/1e3, eff_dis); hold on;
plot3(ess_soc, ones(1,length(ess_soc))*mean(Pdis)/1e3, mean(eff_dis), 'k', 'linewidth', 2);
zlim([0.4 1]);
xlabel('SOC (-)');
ylabel('Discharge Power (kW)');
zlabel('Discharge Efficiency (-)');
view(-120, 30);
defaultratio;
set(gcf, 'pos', [6.6771    5.7708    4.0000    3.0000]);

figure(3); clf;
[C,h] = contour(ess_soc, Pdis/1e3, eff_dis);
caxis([0.5 1]);
clabel(C,h, 'fontsize', 7);
xlabel('SOC (-)');
ylabel('Discharge Power (kW)');
title('Discharge Efficiency (-)');
defaultratio;
set(gcf, 'pos', [10.8750    5.7708    4.0000    3.0000]);


% ================================================================== %%
Pchg = -(1:30)*1e3; % [kW], terminal power, i.e. output power
[SOC,P_CHG] = meshgrid(ess_soc, Pchg);
VOC = interp1(ess_soc, ess_voc, SOC);
R_CHG = interp1(ess_soc, ess_r_chg, SOC);

I_CHG = (VOC - sqrt(VOC.^2 - 4*P_CHG.*R_DIS))./(2*R_DIS);
I_CHG(abs(imag(I_CHG))>0) = nan;
Ptotal = VOC.*I_CHG;
% Ploss = (I_CHG.^2).*R_CHG;

eff_chg = Ptotal./P_CHG;

ave_eff_chg = nanmean(eff_chg(:))

figure(4); clf; 
mesh(ess_soc, Pchg/1e3, eff_chg); hold on;
plot3(ess_soc, ones(1,length(ess_soc))*mean(Pchg)/1e3, mean(eff_chg), 'k', 'linewidth', 2);
zlim([0.4 1]);
xlabel('SOC (-)');
ylabel('Charge Power (kW)');
zlabel('Charge Efficiency (-)');
view(60, 30);
defaultratio;
set(gcf, 'pos', [6.6771    1.7396    4.0000    3.0000]);

figure(5); clf;
[C,h] = contour(ess_soc, Pchg/1e3, eff_chg);
caxis([0.5 1]);
clabel(C,h, 'fontsize', 7);
xlabel('SOC (-)');
ylabel('Charge Power (kW)');
title('Charge Efficiency (-)');
defaultratio;
set(gcf, 'pos', [10.8750    1.7396    4.0000    3.0000]);

