clear all;
close all;

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
load alpha_text;

color_code = rainbow_color(length(alpha_range));
figure(1); clf;
subplot(2,1,1);
hp = plot(1:6, squeeze(J1_collection(1,:,:))/full_potential, 'x-');
for i = 1:length(alpha_range), set(hp(i), 'color', color_code(i,:)); end
ylim([0 1]);
set(gca, 'fontsize', 8);
xlabel('Battery Capacity');
legend(num2str(alpha_range'), 4);
title('J1');

subplot(2,1,2);
score = sum(squeeze(J1_collection(1,:,:)/full_potential));
plot(alpha_range, score, 'x-');
set(gca, 'fontsize', 8);
xlabel('alpha range');
ylabel('total benefit');

% % ==============================
% figure(2); clf;
% subplot(2,1,1);
% plot(1:6, squeeze(J2_collection(1,:,:))/full_potential);
% ylim([0 1]);
% set(gca, 'fontsize', 8);
% xlabel('Battery Capacity');
% legend(num2str(alpha_range'), 4);
% title('J2');
% 
% subplot(2,1,2);
% score = sum(squeeze(J2_collection(1,:,:)/full_potential));
% plot(alpha_range, score, 'x-');
% set(gca, 'fontsize', 8);
% xlabel('alpha range');
% ylabel('total benefit');


%% 0.7:0.05:0.9
load alpha_text_a;

color_code = rainbow_color(length(alpha_range));
figure(101); clf;
subplot(2,1,1);
hp = plot(1:6, squeeze(J1_collection(1,:,:))/full_potential, 'x-');
for i = 1:length(alpha_range), set(hp(i), 'color', color_code(i,:)); end
ylim([0 1]);
set(gca, 'fontsize', 8);
xlabel('Battery Capacity');
legend(num2str(alpha_range'), 4);
title('J1');

subplot(2,1,2);
score = sum(squeeze(J1_collection(1,:,:)/full_potential));
plot(alpha_range, score, 'x-');
set(gca, 'fontsize', 8);
xlabel('alpha range');
ylabel('total benefit');


%% fine
load alpha_text_b;

color_code = rainbow_color(length(alpha_range));
figure(102); clf;
subplot(2,1,1);
hp = plot(1:6, squeeze(J1_collection(1,:,:))/full_potential, 'x-');
for i = 1:length(alpha_range), set(hp(i), 'color', color_code(i,:)); end
ylim([0 1]);
set(gca, 'fontsize', 8);
xlabel('Battery Capacity');
legend(num2str(alpha_range'), 4);
title('J1');

subplot(2,1,2);
score = sum(squeeze(J1_collection(1,:,:)/full_potential));
plot(alpha_range, score, 'x-');
set(gca, 'fontsize', 8);
xlabel('alpha range');
ylabel('total benefit');

