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
figure(1); clf;
plot(1:6, squeeze(J1_collection(1,:,1:6))/full_potential, 'x-');
ylim([0 1]);
xlabel('Battery Capacity');
legend('0','0.2','0.4','0.6','0.8','1', 4);
title('J1');

figure(2); clf;
plot(1:6, squeeze(J2_collection(1,:,1:6))/full_potential);
ylim([0 1]);
xlabel('Battery Capacity');
legend('0','0.2','0.4','0.6','0.8','1', 4);
title('J2');

%%
load alpha_text_finetune;

color_code = rainbow_color(length(alpha_range));
figure(100); clf;
hp = plot(1:6, squeeze(J1_collection(1,:,:))/full_potential, 'x-');
for i = 1:length(alpha_range), set(hp(i), 'color', color_code(i,:)); end
ylim([0 1]);
xlabel('Battery Capacity');
legend(num2str(alpha_range'), 4);
title('J1');

figure(200); clf;
hp = plot(1:6, squeeze(J2_collection(1,:,:))/full_potential);
for i = 1:length(alpha_range), set(hp(i), 'color', color_code(i,:)); end
ylim([0 1]);
xlabel('Battery Capacity');
legend(num2str(alpha_range'), 4);
title('J2');

