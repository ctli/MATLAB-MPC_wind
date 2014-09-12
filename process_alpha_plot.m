clear all;
close all;
clc;
format compact;

data = importdata('SITE_3939_MAIN_FORECASTS.csv', ',');
obs = data.data(:,3);
fcst = data.data(:,4);

scale = 800/max(fcst);
obs = obs*scale;
fcst = fcst*scale;

obs = obs(1:8760)';
fcst = fcst(1:8760)';

c1 = 1;
c2 = 1.03;
c3 = 1;
full_potential = -c1*sum(obs);

%% 0-1, coarse
eta_index = 1; eta = sqrt(0.75);
% eta_index = 2; eta = sqrt(0.95);

bty_size = [0.25 0.5 0.75 1 2 5];

load alpha_text;

color_code = rainbow_color(length(alpha_range));
figure(1); clf;
hp = plot(bty_size, squeeze(J1_collection(eta_index,:,:))/full_potential, 'x-', 'markersize', 3);
for i = 1:length(alpha_range), set(hp(i), 'color', color_code(i,:)); end
ylim([0 1.01]);
set(gca, 'fontsize', 8);
xlim([0.1 5.2]);
set(gca, 'xtick', bty_size, 'xticklabel', {'200 ','400 ','600 ','800 ','1600 ','4000 '});
xlabel('Battery Capacity (MWh)', 'pos', [2.6124   -0.2544   17.3205]);
ylabel('Revenue (Normalized)');
set(gcf, 'unit', 'inch', 'pos', [8.3229    2.1979    3.5000    2.0]);
set(gca, 'units', 'normalized', 'pos', [0.1429    0.2396    0.8065    0.7187]);
set(gca, 'tickdir', 'out');
box off;
rotateticklabel(gca, 90);

my_gridline('y');

legend_string = {'\alpha = 0', '\alpha = 0.1', '\alpha = 0.2', '\alpha = 0.3', '\alpha = 0.4', '\alpha = 0.5', '\alpha = 0.6', '\alpha = 0.7', '\alpha = 0.8', '\alpha = 0.9', '\alpha = 1'};
ax1 = gca;
legend(ax1,flipud(hp(1:5)), fliplr(legend_string(1:5)));
set(legend, 'units', 'pixels');
set(legend, 'pos', [187.0000   52   53.8333   65.8333]);
set(legend, 'box', 'off');

ah2=axes('position',get(gca,'position'), 'visible','off');
legend(ah2,flipud(hp(6:end)),  fliplr(legend_string(6:end)));
set(legend,'units','pixels', 'fontsize', 8);
set(legend, 'pos', [259.8333   51   53.8333   80.8333]);
set(legend, 'box', 'off');

% set(gcf, 'color', 'w');
% export_fig alpha -painters;


%%
figure(2); clf; hold on;
score = mean(squeeze(J1_collection(eta_index,:,:)/full_potential));
plot(alpha_range, score, 'k-');
for i = 1:11
    plot(alpha_range(i), score(i), 'x', 'color', color_code(i,:));
end;
ylim([0.3 0.503]);
set(gca, 'fontsize', 8);
xlabel('\alpha', 'pos', [0.4926    0.269   17.3205]);
ylabel('Averaged Revenue (Normalized)');
set(gcf, 'unit', 'inch', 'pos', [12.0208    2.1979    3.5000    1.85]);
set(gca, 'units', 'pixels', 'pos', [49.0144   32.5238  270.9840  134.4762]);
box off;
set(gca, 'tickdir', 'out');
my_gridline('y');
% export_fig alpha_range -painters;


