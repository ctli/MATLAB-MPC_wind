clear all
close all
clc
format compact

r = {'BJ','TJ','HE','SX','SD','NM','LN','JL','HL','SH','JS','ZJ','AH','FJ','HA','HB','HN','JX','SC','CQ','SN','GS','QH','NX','XJ','GD','GX','YN','GZ','HI'};
g = {'AGR', 'COL', 'GAS', 'CRU', 'OIL', 'ELE', 'EIS', 'MAN', 'OMN', 'CON', 'TRN', 'WTR', 'SER', 'c', 'OTF'};
yr = [2007,2010,2015,2020,2025,2030];


% %% Urban emission breakdown (combustion emission vs. process emission)
% [urban_c, urban_c_id] = getgdx('result_urban_exo.gdx', 'urban_c'); % [urb]x[rs]x[g]x[t] = [9]x[30]x[8]x[6]
% [urban_o, urban_o_id] = getgdx('result_urban_exo.gdx', 'urban_o'); % [urb]x[rs]x[g]x[t] = [9]x[30]x[10]x[6]
% 
% fc1 = rainbow_color(15);
% fc2 = flipud(rainbow_color(16));
% fc2(1,:) = [];
% fc = [fc1;fc2];
% 
% for i = [1, 2, 5]%:length(urban_id{1})
% urb_c = squeeze(urban_c(strcmp(urban_c_id{1}(i), urban_c_id{1}),:,:,:)); % [30]x[8]x[6]
% urb_c_2007 = urb_c(:,:,1); % [30]x[8]; AGR/OMN/EIS/MAN/ELE/TRN/SER/c
% 
% urb_o = squeeze(urban_o(strcmp(urban_o_id{1}(i), urban_o_id{1}),:,:,:)); % [30]x[8]x[6]
% urb_o = urb_o(:,:,1); % [30]x[8]; AGR/OMN/EIS/MAN/ELE/TRN/SER/c
% 
% temp1 = zeros(30, length(g));
% temp2 = zeros(30, length(g));
% 
% [~,id] = ismember(urban_c_id{3}, g);
% temp1(:,id) = urb_c_2007;
% 
% [~,id] = ismember(urban_o_id{3}, g);
% temp2(:,id) = urb_o;
% 
% figure(i); clf;
% hb0 = bar((1:length(g))-0.16, temp1', 0.3, 'stacked', 'edge', 'k'); hold on;
% for j = 1:length(hb0), set(hb0(j), 'facec', fc(j,:), 'edge', 'k'); end
% hb1 = bar((1:length(g))-0.16, sum(temp1), 0.3); hold on;
% set(hb1, 'facec', 'none', 'edge', 'r', 'linewidth', 1);
% 
% hb0 = bar((1:length(g))+0.16, temp2', 0.3, 'stacked', 'edge', 'k');
% for j = 1:length(hb0), set(hb0(j), 'facec', fc(j,:), 'edge', 'k'); end
% hb2 = bar((1:length(g))+0.16, sum(temp2), 0.3, 'stacked');
% set(hb2, 'facec', 'none', 'edge', [0 0.6 0], 'linewidth', 1);
% 
% xlim([0 16]);
% set(gca, 'xtick', 1:length(g), 'xticklabel', g, 'fontsize', 9);
% rotateticklabel(gca, 90);
% ylabel([char(cellstr(urban_c_id{1}(i))), ' in 2007 (Tg)'], 'fontsize', 12, 'fontweight', 'bold');
% legend([hb1, hb2], 'combustion', 'process', 'location', 'north');
% 
% if i == 2
%     for rn = 1:length(r)
%     text(15.5, sum(temp2(1:rn-1,end))+temp2(rn,end)/2, r(rn), 'fontsize', 6);
%     end
% end
% 
% set(gcf, 'units', 'inch', 'pos', [1 1 3.5, 5.85]);
% set(gca, 'pos', [0.1667    0.1100    0.7562    0.8150]);
% 
% export_name = [char(cellstr(urban_c_id{1}(i))), '_sectorial'];
% % export_fig(export_name);
% end


%% Urban emission breakdown (combustion emission vs. process emission)
[urban_c, urban_c_id] = getgdx('result_urban_exo.gdx', 'urban_c'); % [urb]x[rs]x[g]x[t] = [9]x[30]x[8]x[6]
[urban_o, urban_o_id] = getgdx('result_urban_exo.gdx', 'urban_o'); % [urb]x[rs]x[g]x[t] = [9]x[30]x[11]x[6]
for i = 1:length(urban_c_id{1})
urb_c = squeeze(urban_c(strcmp(urban_c_id{1}(i), urban_c_id{1}),:,:,:)); % [30]x[8]x[6]
urb_c_g = squeeze(sum(urb_c,2)); % [30]x[6]
urb_c_n = sum(urb_c_g); % [1]x[6]

urb_o = squeeze(urban_o(strcmp(urban_o_id{1}(i), urban_o_id{1}),:,:,:)); % [30]x[8]x[6]
urb_o_g = squeeze(sum(urb_o(:,1:end-1,:),2)); % [30]x[6]
urb_o_n = sum(urb_o_g); % [1]x[6]
urb_o_b = squeeze(sum(urb_o(:,end,:))); % biomass burning

figure(400+i); clf; hold on; box on;
hb0 = bar(yr-0.65, [urb_c_n; urb_o_n; urb_o_b']', 0.25, 'stacked', 'edge', 'none');
% set(hb0(1), 'facec', [0.8 0 0]); %combustion
% set(hb0(2), 'facec', [1 0.65 0]); %process
% hb = bar(yr-0.5, [urb_c_n; urb_o_n; urb_o_b']', 0.3, 'stacked', 'edge', [0.5 0.5 0.5], 'linewidth', 2);
% set(hb(1), 'facec', [0.8 0 0]); %combustion
% set(hb(2), 'facec', [1 0.65 0]); %process
% set(hb(3), 'facec', [1 1 0]); %biomass burning

hb1 = bar(yr-0.65, sum([urb_c_n; urb_o_n; urb_o_b']), 0.25, 'stacked', 'facec', 'none', 'edge', [0.5 0.5 0.5], 'linewidth', 1);

colormap copper
set(hb0(1), 'facec', [0.5 0 0]); %biomass burning

set(gca, 'fontsize', 10);
set(gca, 'layer', 'top');
set(gcf, 'unit', 'inch', 'pos', [0.25+4.15*(i-1)    5.25    4.0000    3.0000]);
set(gca, 'pos', [0.1458    0.1100    0.7722    0.8150]);
set(gca, 'layer', 'bottom');
ylabel([char(cellstr(urban_c_id{1}(i))), ' (Tg)'], 'fontsize', 12, 'fontweight', 'bold');
xlim([2005 2031.5]);
legend(fliplr(hb0), 'Biomass Burning & Other', 'Industrial Process', 'Fossil Fuel Combustion', 2);
set(legend, 'fontsize', 8);

if i == 1, ylim([0 3]); end
if i == 2, ylim([0 500]); end
if i == 3, ylim([0 30]); end
if i == 4, ylim([0 80]); end
if i == 5, ylim([0 6]); end
if i == 6, ylim([0 55]); end
if i == 7, ylim([0 40]); end
if i == 8, ylim([0 70]); end
if i == 9, ylim([0 65]); end
end

%%
[urban_c, urban_c_id] = getgdx('result_egyint_n.gdx', 'urban_c'); % [urb]x[rs]x[g]x[t] = [9]x[30]x[8]x[6]
[urban_o, urban_o_id] = getgdx('result_egyint_n.gdx', 'urban_o'); % [urb]x[rs]x[g]x[t] = [9]x[30]x[10]x[6]
for i = 1:length(urban_c_id{1})
urb_c = squeeze(urban_c(strcmp(urban_c_id{1}(i), urban_c_id{1}),:,:,:)); % [30]x[8]x[6]
urb_c_g = squeeze(sum(urb_c,2)); % [30]x[6]
urb_c_n = sum(urb_c_g); % [1]x[6]

urb_o = squeeze(urban_o(strcmp(urban_o_id{1}(i), urban_o_id{1}),:,:,:)); % [30]x[8]x[6]
urb_o_g = squeeze(sum(urb_o(:,1:end-1,:),2)); % [30]x[6]
urb_o_n = sum(urb_o_g); % [1]x[6]
urb_o_b = squeeze(sum(urb_o(:,end,:))); % biomass burning

figure(400+i); 
hb0 = bar(yr+0, [urb_c_n; urb_o_n; urb_o_b']', 0.25, 'stacked', 'edge', 'none');
set(hb0(1), 'facec', [0.5 0 0]); %biomass burning
hb1 = bar(yr+0, sum([urb_c_n; urb_o_n; urb_o_b']), 0.25, 'stacked', 'facec', 'none', 'edge', [0 0.6 0], 'linewidth', 1);

file_name = ['fig_', char(cellstr(urban_c_id{1}(i)))];
% export_fig(file_name);
end

%%
[urban_c, urban_c_id] = getgdx('result_ccap_r.gdx', 'urban_c'); % [urb]x[rs]x[g]x[t] = [9]x[30]x[8]x[6]
[urban_o, urban_o_id] = getgdx('result_ccap_r.gdx', 'urban_o'); % [urb]x[rs]x[g]x[t] = [9]x[30]x[10]x[6]
for i = 1:length(urban_c_id{1})
urb_c = squeeze(urban_c(strcmp(urban_c_id{1}(i), urban_c_id{1}),:,:,:)); % [30]x[8]x[6]
urb_c_g = squeeze(sum(urb_c,2)); % [30]x[6]
urb_c_n = sum(urb_c_g); % [1]x[6]

urb_o = squeeze(urban_o(strcmp(urban_o_id{1}(i), urban_o_id{1}),:,:,:)); % [30]x[8]x[6]
urb_o_g = squeeze(sum(urb_o(:,1:end-1,:),2)); % [30]x[6]
urb_o_n = sum(urb_o_g); % [1]x[6]
urb_o_b = squeeze(sum(urb_o(:,end,:))); % biomass burning

figure(400+i); 
hb0 = bar(yr+0.65, [urb_c_n; urb_o_n; urb_o_b']', 0.25, 'stacked', 'edge', 'none');
set(hb0(1), 'facec', [0.5 0 0]); %biomass burning
hb1 = bar(yr+0.65, sum([urb_c_n; urb_o_n; urb_o_b']), 0.25, 'stacked', 'facec', 'none', 'edge', 'b', 'linewidth', 1);

if i == 1
    text(2002.1, 3.15, 'Gray: No Policy; ', 'fontsize', 8, 'color', [0.5 0.5 0.5]);
    text(2009.3, 3.15, 'Green: Ntnl Energy Intensity Target; ', 'fontsize', 8, 'color', [0 0.6 0]);
    text(2025, 3.15, 'Blue: Rgnl CO2 Cap', 'fontsize', 8, 'color', 'b');
end

file_name = ['fig_', char(cellstr(urban_c_id{1}(i)))];
export_fig(file_name, '-painters');
end



% %% compare AGR & private consumption outputs
% % this helps to understand if it is more appropriate to link biomass burning to AGR, instead of c
% [Y, Y_id] = getgdx('merge_urb_exo.gdx', 'Y'); % [t]x[g]x[rs] = [6]x[18]x[34]
% YY = Y(:,1:16,1:30);
% YY_n = sum(YY,3);
% YY_n_AGR = YY_n(:,strcmp(Y_id{2}, 'AGR'));
% YY_n_c   = YY_n(:,strcmp(Y_id{2}, 'c'));
% 
% figure(4); clf;
% plot(yr, YY_n_AGR, 'x-', yr, YY_n_c, 'x-');
% legend('AGR', 'c', 2);
% ylabel('Sectorial Output (Normalized)');
% 
% [vom, vom_id] = getgdx('result_urban_exo.gdx', 'vom'); % [g]x[rs] = [16]x[34]
% value_AGR_n = zeros(1,length(yr));
% for yn = 1:length(yr)
% value_AGR = squeeze(YY(yn,strcmp(Y_id{2}, 'AGR'),1:30)).*vom(strcmp(vom_id{1}, 'AGR'),1:30)';
% value_AGR_n(yn) = sum(value_AGR(:));
% end
% value_c_n = zeros(1,length(yr));
% for yn = 1:length(yr)
% value_c = squeeze(YY(yn,strcmp(Y_id{2}, 'c'),1:30)).*vom(strcmp(vom_id{1}, 'c'),1:30)';
% value_c_n(yn) = sum(value_c(:));
% end
% 
% figure(5); clf;
% plot(yr, value_AGR_n, 'x-', yr, value_c_n, 'x-');
% ylabel('Sectorial Output (Billion USD in 2007)');
% legend('AGR', 'c', 2);

